/*
 * 細胞内分子夾雑シミュレーション（全画面モード）
 *
 * 大腸菌の細胞質から一辺 L nm の立方体を切り出し、その中の高分子を
 * 実際の細胞内濃度のまま配置してブラウン動力学で動かします。
 *
 * 物理:
 *   過減衰ランジュバン方程式（慣性項を無視した拡散支配の運動）
 *     dx = (F/γ)·dt + sqrt(2·D·dt)·ξ,   γ = kT/D  （アインシュタインの関係式）
 *   排除体積は重なり量に比例する斥力（ソフトコア）で表現します。
 *   境界は周期境界条件。平均二乗変位(MSD)は折り返しを展開した座標で測ります。
 *
 * 単位系: 長さ nm, 時間 µs, エネルギー kT（= 1）
 */
(function(){
  'use strict';

  // --- 大腸菌細胞質の組成（1 µm³ あたりのコピー数）-------------------------
  // 半径は流体力学的半径のおおよその値。
  const SPECIES = [
    { key:'ribosome', label:'リボソーム (70S)',   r:10.5, perUm3:20000,   color:'#c084fc' },
    { key:'mrna',     label:'mRNA',               r:8.0,  perUm3:1400,    color:'#fbbf24' },
    { key:'rnap',     label:'RNAポリメラーゼ',     r:5.0,  perUm3:2000,    color:'#60a5fa' },
    { key:'protein',  label:'その他のタンパク質',   r:2.5,  perUm3:2400000, color:'#8794a8' }
  ];
  // 回路から発現した分子はここに追加されます（app.js が渡します）。
  const D_REF_R = 2.4, D_REF = 8.0;      // GFP: r≈2.4nm, D≈8 µm²/s = 8 nm²/µs
  const K_REP = 16;                       // ソフトコア斥力のばね定数 [kT/nm²]
  const SKIN = 2.0;                       // Verlet リストのスキン幅 [nm]
  const DT = 0.005;                       // 時間刻み [µs]

  const diffusion = r => D_REF * (D_REF_R / r);   // ストークス・アインシュタイン D ∝ 1/r

  // =====================================================================
  // シミュレーション本体
  // =====================================================================
  function Sim(){
    this.L = 150;
    this.N = 0;
    this.pos = null; this.unwrapped = null; this.origin = null;
    this.radius = null; this.diff = null; this.spec = null;
    this.speciesList = [];
    this.time = 0;
    this.crowding = true;
    this.pairs = null; this.pairCount = 0; this.stepsSinceList = 1e9;
    this.maxR = 1;
  }

  Sim.prototype.build = function(L, extraSpecies){
    this.L = L;
    const vol = Math.pow(L / 1000, 3);            // µm³
    const list = SPECIES.concat(extraSpecies || []);
    const counts = list.map(s => Math.max(0, Math.round(s.perUm3 * vol)));
    const N = counts.reduce((a, b) => a + b, 0);

    this.speciesList = list.map((s, i) => Object.assign({}, s, { count: counts[i], index: i }));
    this.N = N;
    this.pos = new Float32Array(N * 3);
    this.unwrapped = new Float32Array(N * 3);
    this.radius = new Float32Array(N);
    this.diff = new Float32Array(N);
    this.spec = new Uint8Array(N);

    let k = 0;
    this.speciesList.forEach(s => {
      const D = diffusion(s.r);
      for (let n = 0; n < s.count; n++, k++) {
        this.pos[k*3]   = Math.random() * L;
        this.pos[k*3+1] = Math.random() * L;
        this.pos[k*3+2] = Math.random() * L;
        this.radius[k] = s.r;
        this.diff[k] = D;
        this.spec[k] = s.index;
      }
    });
    this.unwrapped.set(this.pos);
    this.origin = new Float32Array(this.pos);
    this.maxR = Math.max.apply(null, this.speciesList.filter(s => s.count > 0).map(s => s.r).concat([1]));
    this.time = 0;
    this.stepsSinceList = 1e9;
    this.pairs = new Int32Array(Math.max(1024, N * 48));
    // 過密配置の重なりを先に緩和しておく（ペアリストは数ステップおきに更新）
    for (let i = 0; i < 90; i++) {
      if (i % 8 === 0) this.buildPairs();
      this.relax();
    }
    this.unwrapped.set(this.pos);
    this.origin.set(this.pos);
    return this;
  };

  Sim.prototype.volumeFraction = function(){
    let v = 0;
    for (let i = 0; i < this.N; i++) {
      const r = this.radius[i];
      v += 4 / 3 * Math.PI * r * r * r;
    }
    return v / (this.L * this.L * this.L);
  };

  // --- 近接ペアリスト（セルリスト法 + Verlet スキン）----------------------
  //
  // 粒子の大きさが4倍以上違う（リボソーム 10.5nm / タンパク質 2.5nm）ため、
  // 全粒子を最大半径で探索すると数の多い小粒子が無駄に広い範囲を走査します。
  // そこで各粒子の探索半径を 2·r_i + SKIN とし、ペアは「大きい方の粒子」が
  // 登録します。r_i ≥ r_j なら r_i+r_j+SKIN ≤ 2·r_i+SKIN なので取りこぼしません。
  Sim.prototype.buildPairs = function(){
    const L = this.L, N = this.N, pos = this.pos, radius = this.radius;
    const cs0 = 7;                                     // 小粒子スケールのセル幅 [nm]
    const nc = Math.max(1, Math.floor(L / cs0));
    const cs = L / nc;
    const half = L / 2;

    let head = this._head;
    if (!head || head.length !== nc*nc*nc) head = this._head = new Int32Array(nc*nc*nc);
    head.fill(-1);
    const next = this._next || (this._next = new Int32Array(N));
    const ci = this._ci || (this._ci = new Int32Array(N));
    const cj = this._cj || (this._cj = new Int32Array(N));
    const ck = this._ck || (this._ck = new Int32Array(N));

    for (let i = 0; i < N; i++) {
      const a = Math.min(nc - 1, Math.floor(pos[i*3]   / cs));
      const b = Math.min(nc - 1, Math.floor(pos[i*3+1] / cs));
      const c = Math.min(nc - 1, Math.floor(pos[i*3+2] / cs));
      ci[i] = a; cj[i] = b; ck[i] = c;
      const h = (a * nc + b) * nc + c;
      next[i] = head[h]; head[h] = i;
    }

    let p = 0;
    const cap = this.pairs.length - 2;
    const maxSpan = Math.max(0, (nc - 1) >> 1);
    for (let i = 0; i < N; i++) {
      const ri = radius[i];
      const span = Math.min(maxSpan, Math.ceil((2 * ri + SKIN) / cs));
      const xi = pos[i*3], yi = pos[i*3+1], zi = pos[i*3+2];
      for (let da = -span; da <= span; da++) {
        let a = (ci[i] + da) % nc; if (a < 0) a += nc;
        for (let db = -span; db <= span; db++) {
          let b = (cj[i] + db) % nc; if (b < 0) b += nc;
          for (let dc = -span; dc <= span; dc++) {
            let c = (ck[i] + dc) % nc; if (c < 0) c += nc;
            let j = head[(a * nc + b) * nc + c];
            while (j !== -1) {
              const rj = radius[j];
              // 大きい方が登録する（同径なら添字の小さい方）
              if (rj < ri || (rj === ri && i < j)) {
                let dx = pos[j*3] - xi, dy = pos[j*3+1] - yi, dz = pos[j*3+2] - zi;
                if (dx >  half) dx -= L; else if (dx < -half) dx += L;
                if (dy >  half) dy -= L; else if (dy < -half) dy += L;
                if (dz >  half) dz -= L; else if (dz < -half) dz += L;
                const cut = ri + rj + SKIN;
                if (dx*dx + dy*dy + dz*dz < cut*cut && p < cap) {
                  this.pairs[p++] = i; this.pairs[p++] = j;
                }
              }
              j = next[j];
            }
          }
        }
      }
    }
    this.pairCount = p >> 1;
    this.stepsSinceList = 0;
  };

  // 初期配置の重なりを取るための位置ずらし（熱ゆらぎ無し）
  Sim.prototype.relax = function(){
    const L = this.L, pos = this.pos, radius = this.radius, pairs = this.pairs;
    for (let n = 0; n < this.pairCount; n++) {
      const i = pairs[n*2], j = pairs[n*2+1];
      let dx = pos[j*3] - pos[i*3], dy = pos[j*3+1] - pos[i*3+1], dz = pos[j*3+2] - pos[i*3+2];
      if (dx >  L/2) dx -= L; else if (dx < -L/2) dx += L;
      if (dy >  L/2) dy -= L; else if (dy < -L/2) dy += L;
      if (dz >  L/2) dz -= L; else if (dz < -L/2) dz += L;
      const d2 = dx*dx + dy*dy + dz*dz, s = radius[i] + radius[j];
      if (d2 < s*s && d2 > 1e-9) {
        const d = Math.sqrt(d2), push = (s - d) * 0.25 / d;
        const px = dx*push, py = dy*push, pz = dz*push;
        pos[i*3] -= px; pos[i*3+1] -= py; pos[i*3+2] -= pz;
        pos[j*3] += px; pos[j*3+1] += py; pos[j*3+2] += pz;
      }
    }
    this.wrap();
  };

  Sim.prototype.wrap = function(){
    const L = this.L, pos = this.pos;
    for (let i = 0; i < this.N * 3; i++) {
      if (pos[i] < 0) pos[i] += L; else if (pos[i] >= L) pos[i] -= L;
    }
  };

  let gaussSpare = null;
  function gauss(){
    if (gaussSpare !== null) { const v = gaussSpare; gaussSpare = null; return v; }
    let u, v, s;
    do { u = Math.random()*2-1; v = Math.random()*2-1; s = u*u + v*v; } while (s >= 1 || s === 0);
    const f = Math.sqrt(-2 * Math.log(s) / s);
    gaussSpare = v * f;
    return u * f;
  }

  Sim.prototype.step = function(){
    const L = this.L, N = this.N, pos = this.pos, un = this.unwrapped;
    const radius = this.radius, diff = this.diff;

    if (this.stepsSinceList > 8) this.buildPairs();

    // 排除体積による力
    const fx = this._fx || (this._fx = new Float32Array(N));
    const fy = this._fy || (this._fy = new Float32Array(N));
    const fz = this._fz || (this._fz = new Float32Array(N));
    fx.fill(0); fy.fill(0); fz.fill(0);

    if (this.crowding) {
      const pairs = this.pairs;
      for (let n = 0; n < this.pairCount; n++) {
        const i = pairs[n*2], j = pairs[n*2+1];
        let dx = pos[j*3] - pos[i*3], dy = pos[j*3+1] - pos[i*3+1], dz = pos[j*3+2] - pos[i*3+2];
        if (dx >  L/2) dx -= L; else if (dx < -L/2) dx += L;
        if (dy >  L/2) dy -= L; else if (dy < -L/2) dy += L;
        if (dz >  L/2) dz -= L; else if (dz < -L/2) dz += L;
        const d2 = dx*dx + dy*dy + dz*dz, s = radius[i] + radius[j];
        if (d2 < s*s && d2 > 1e-9) {
          const d = Math.sqrt(d2), mag = K_REP * (s - d) / d;
          const ax = dx*mag, ay = dy*mag, az = dz*mag;
          fx[i] -= ax; fy[i] -= ay; fz[i] -= az;
          fx[j] += ax; fy[j] += ay; fz[j] += az;
        }
      }
    }

    // 過減衰ランジュバン: dx = D·F/kT·dt + sqrt(2·D·dt)·ξ   (kT = 1)
    for (let i = 0; i < N; i++) {
      const D = diff[i], mob = D * DT, amp = Math.sqrt(2 * D * DT);
      const dx = fx[i]*mob + gauss()*amp;
      const dy = fy[i]*mob + gauss()*amp;
      const dz = fz[i]*mob + gauss()*amp;
      un[i*3] += dx; un[i*3+1] += dy; un[i*3+2] += dz;
      let x = pos[i*3] + dx, y = pos[i*3+1] + dy, z = pos[i*3+2] + dz;
      if (x < 0) x += L; else if (x >= L) x -= L;
      if (y < 0) y += L; else if (y >= L) y -= L;
      if (z < 0) z += L; else if (z >= L) z -= L;
      pos[i*3] = x; pos[i*3+1] = y; pos[i*3+2] = z;
    }
    this.stepsSinceList++;
    this.time += DT;
  };

  // 種ごとの平均二乗変位から見かけの拡散係数を出す: MSD = 6·D·t
  Sim.prototype.msdOf = function(speciesIndex){
    const un = this.unwrapped, or = this.origin;
    let sum = 0, n = 0;
    for (let i = 0; i < this.N; i++) {
      if (this.spec[i] !== speciesIndex) continue;
      const dx = un[i*3] - or[i*3], dy = un[i*3+1] - or[i*3+1], dz = un[i*3+2] - or[i*3+2];
      sum += dx*dx + dy*dy + dz*dz; n++;
    }
    return n ? sum / n : 0;
  };

  Sim.prototype.resetMSD = function(){
    this.origin.set(this.unwrapped);
    this.time = 0;
  };

  // =====================================================================
  // 全画面ビュー
  // =====================================================================
  let overlay = null, sim = null, renderer = null, scene = null, camera = null;
  let meshes = [], raf = null, running = false;
  let orbit = { theta: 0.8, phi: 1.1, dist: 380 };
  let dragging = false, lastPtr = null;
  let extraSpecies = [], speedSubsteps = 3, fpsEMA = 60, lastFrame = 0;

  function el(tag, cls, html){
    const e = document.createElement(tag);
    if (cls) e.className = cls;
    if (html != null) e.innerHTML = html;
    return e;
  }

  function buildScene(){
    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(42, 1, 1, 4000);
    if (window.MolLook) { window.MolLook.setupScene(scene, 0.002); return; }
    scene.background = new THREE.Color(0x070d12);
    scene.add(new THREE.AmbientLight(0xffffff, 0.62));
    const d1 = new THREE.DirectionalLight(0xffffff, 0.85); d1.position.set(1, 1, 1); scene.add(d1);
    const d2 = new THREE.DirectionalLight(0x7fb6ff, 0.35); d2.position.set(-1, -0.6, -0.8); scene.add(d2);
  }

  function rebuildMeshes(){
    meshes.forEach(m => { scene.remove(m.mesh); m.mesh.geometry.dispose(); m.mesh.material.dispose(); });
    meshes = [];
    const L = sim.L;
    sim.speciesList.forEach(s => {
      if (!s.count) return;
      // 小さい粒子ほど数が多いのでポリゴンを落とす
      const detail = s.r > 7 ? 2 : (s.r > 4 ? 1 : 0);
      const geo = new THREE.IcosahedronGeometry(s.r, detail);
      const mat = window.MolLook
        ? window.MolLook.material(s.color, { roughness: 0.6, rimStr: s.key === 'protein' ? 0.35 : 0.6,
            transparent: s.key === 'protein', opacity: s.key === 'protein' ? 0.85 : 1 })
        : new THREE.MeshStandardMaterial({
            color: new THREE.Color(s.color), roughness: 0.5, metalness: 0.05,
            transparent: s.key === 'protein', opacity: s.key === 'protein' ? 0.85 : 1
          });
      const mesh = new THREE.InstancedMesh(geo, mat, s.count);
      mesh.frustumCulled = false;
      mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
      // 単位行列を並べておき、以後は平行移動成分だけ書き換える
      const arr = mesh.instanceMatrix.array;
      for (let n = 0; n < s.count; n++) {
        const o = n * 16;
        arr[o] = 1; arr[o+5] = 1; arr[o+10] = 1; arr[o+15] = 1;
      }
      scene.add(mesh);
      meshes.push({ mesh, species: s, visible: true });
    });

    const box = new THREE.LineSegments(
      new THREE.EdgesGeometry(new THREE.BoxGeometry(L, L, L)),
      new THREE.LineBasicMaterial({ color: 0x2f4f66 })
    );
    box.name = '__box';
    const prev = scene.getObjectByName('__box');
    if (prev) scene.remove(prev);
    scene.add(box);
  }

  function updateInstances(){
    const L = sim.L, half = L / 2, pos = sim.pos, spec = sim.spec;
    const cursor = new Int32Array(meshes.length);
    const byIndex = {};
    meshes.forEach((m, i) => { byIndex[m.species.index] = i; });
    for (let i = 0; i < sim.N; i++) {
      const mi = byIndex[spec[i]];
      if (mi === undefined) continue;
      const m = meshes[mi], n = cursor[mi]++;
      const o = n * 16, arr = m.mesh.instanceMatrix.array;
      arr[o+12] = pos[i*3]   - half;
      arr[o+13] = pos[i*3+1] - half;
      arr[o+14] = pos[i*3+2] - half;
    }
    meshes.forEach(m => { m.mesh.instanceMatrix.needsUpdate = true; });
  }

  function frame(ts){
    if (!running) return;
    for (let s = 0; s < speedSubsteps; s++) sim.step();
    updateInstances();

    const cx = orbit.dist * Math.sin(orbit.phi) * Math.cos(orbit.theta);
    const cy = orbit.dist * Math.cos(orbit.phi);
    const cz = orbit.dist * Math.sin(orbit.phi) * Math.sin(orbit.theta);
    camera.position.set(cx, cy, cz);
    camera.lookAt(0, 0, 0);
    renderer.render(scene, camera);

    if (lastFrame) fpsEMA = fpsEMA * 0.9 + (1000 / Math.max(1, ts - lastFrame)) * 0.1;
    lastFrame = ts;
    updateStats();
    raf = requestAnimationFrame(frame);
  }

  // 拡散係数は短時間側では自由拡散に漸近するため、長時間極限に入るまでは
  // 値を出さずに測定中と表示します（MSD = 6·D·t が線形になる領域で評価）。
  const MSD_SETTLE = 2.0;   // [µs]

  function updateStats(){
    const box = overlay.querySelector('#cwStats');
    if (!box) return;
    const proteinIdx = sim.speciesList.findIndex(s => s.key === 'protein');
    const dApp = sim.time > 0 ? sim.msdOf(proteinIdx) / (6 * sim.time) : 0;
    const dFree = diffusion(2.5);
    const settled = sim.time >= MSD_SETTLE;
    const pending = '<span class="cw-pending">測定中… ' +
      Math.max(0, MSD_SETTLE - sim.time).toFixed(1) + ' µs</span>';
    box.innerHTML =
      row('分子数', sim.N.toLocaleString() + ' 個') +
      row('体積占有率', (sim.volumeFraction() * 100).toFixed(1) + ' %') +
      row('経過時間', sim.time.toFixed(2) + ' µs') +
      row('見かけの D', settled ? dApp.toFixed(2) + ' nm²/µs' : pending) +
      row('自由拡散比', settled ? (dApp / dFree * 100).toFixed(0) + ' %' : pending) +
      row('描画', fpsEMA.toFixed(0) + ' fps');
  }
  const row = (k, v) => '<div class="cw-stat"><span>' + k + '</span><b>' + v + '</b></div>';

  function buildLegend(){
    const host = overlay.querySelector('#cwLegend');
    host.innerHTML = '';
    sim.speciesList.forEach(s => {
      if (!s.count) return;
      const item = el('label', 'cw-legend-item');
      item.innerHTML =
        '<input type="checkbox" checked data-key="' + s.key + '">' +
        '<i style="background:' + s.color + '"></i>' +
        '<span class="cw-legend-name">' + s.label + '</span>' +
        '<span class="cw-legend-count">' + s.count.toLocaleString() + '</span>' +
        '<span class="cw-legend-r">r=' + s.r + 'nm</span>';
      host.appendChild(item);
      item.querySelector('input').addEventListener('change', e => {
        const m = meshes.find(x => x.species.key === s.key);
        if (m) m.mesh.visible = e.target.checked;
      });
    });
  }

  function resize(){
    if (!renderer) return;
    const host = overlay.querySelector('#cwCanvas');
    const w = host.clientWidth, h = host.clientHeight;
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 1.6));
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  }

  function rebuild(L){
    sim.build(L, extraSpecies);
    orbit.dist = L * 2.5;
    // 箱の手前と奥で霧の濃さが変わり、奥行きが分かる程度に
    if (scene.fog) scene.fog.density = 0.6 / (L * 2.5);
    rebuildMeshes();
    buildLegend();
    updateInstances();
    drawInset();
  }

  function open(opts){
    if (typeof THREE === 'undefined') {
      alert('3D描画ライブラリ(Three.js)を読み込めませんでした。インターネット接続を確認してください。');
      return;
    }
    extraSpecies = (opts && opts.extraSpecies) || [];
    if (!overlay) overlay = buildOverlay();
    overlay.style.display = 'flex';
    document.body.style.overflow = 'hidden';

    if (!renderer) {
      buildScene();
      renderer = new THREE.WebGLRenderer({ antialias: true, powerPreference: 'high-performance' });
      if (window.MolLook) window.MolLook.setupRenderer(renderer);
      overlay.querySelector('#cwCanvas').appendChild(renderer.domElement);
      attachOrbit(renderer.domElement);
      window.addEventListener('resize', resize);
    }
    sim = new Sim();
    const L = Number(overlay.querySelector('#cwSize').value);
    rebuild(L);
    resize();
    running = true; lastFrame = 0;
    raf = requestAnimationFrame(frame);
  }

  function close(){
    running = false;
    if (raf) cancelAnimationFrame(raf);
    raf = null;
    if (overlay) overlay.style.display = 'none';
    document.body.style.overflow = '';
  }

  function attachOrbit(dom){
    dom.addEventListener('pointerdown', e => { dragging = true; lastPtr = { x:e.clientX, y:e.clientY }; });
    window.addEventListener('pointermove', e => {
      if (!dragging) return;
      orbit.theta -= (e.clientX - lastPtr.x) * 0.006;
      orbit.phi = Math.max(0.12, Math.min(Math.PI - 0.12, orbit.phi - (e.clientY - lastPtr.y) * 0.006));
      lastPtr = { x:e.clientX, y:e.clientY };
    });
    window.addEventListener('pointerup', () => { dragging = false; });
    dom.addEventListener('wheel', e => {
      e.preventDefault();
      orbit.dist = Math.max(sim.L * 0.6, Math.min(sim.L * 6, orbit.dist + e.deltaY * 0.4));
    }, { passive:false });
  }

  function buildOverlay(){
    const o = el('div', 'cw-overlay');
    o.innerHTML =
      '<div class="cw-head">' +
        '<div class="cw-title">細胞内分子夾雑シミュレーション' +
          '<small>大腸菌の細胞質から切り出した立方体を、実際の細胞内濃度のままブラウン動力学で計算しています</small>' +
        '</div>' +
        '<button class="cw-close" id="cwClose">閉じる ✕</button>' +
      '</div>' +
      '<div class="cw-body">' +
        '<div class="cw-canvas" id="cwCanvas"><div class="cw-hint">ドラッグで回転 / ホイールでズーム</div></div>' +
        '<aside class="cw-side">' +
          '<div class="cw-panel">' +
            '<h4>サンプリング体積</h4>' +
            '<div class="cw-inset" id="cwInset"></div>' +
            '<div class="cw-row"><label>一辺</label>' +
              '<input type="range" id="cwSize" min="80" max="260" step="10" value="150">' +
              '<span id="cwSizeVal">150 nm</span></div>' +
            '<div class="cw-note" id="cwCountNote"></div>' +
          '</div>' +
          '<div class="cw-panel">' +
            '<h4>計測値</h4>' +
            '<div id="cwStats"></div>' +
            '<button class="cw-btn" id="cwResetMsd">MSDを測り直す</button>' +
          '</div>' +
          '<div class="cw-panel">' +
            '<h4>操作</h4>' +
            '<div class="cw-row"><label>速さ</label>' +
              '<input type="range" id="cwSpeed" min="1" max="10" step="1" value="3">' +
              '<span id="cwSpeedVal">×3</span></div>' +
            '<label class="cw-check"><input type="checkbox" id="cwCrowd" checked> 排除体積（分子夾雑）を有効にする</label>' +
            '<div class="cw-note">オフにすると分子が互いをすり抜けます（自由拡散）。切り替えるとMSDの測定は自動でやり直されます。<b>排除体積だけ</b>でも拡散は3割ほど遅くなりますが、実際の細胞ではさらに非特異的な結合などが加わり、もっと遅くなることが知られています。</div>' +
          '</div>' +
          '<div class="cw-panel">' +
            '<h4>分子種</h4>' +
            '<div id="cwLegend"></div>' +
          '</div>' +
        '</aside>' +
      '</div>';
    document.body.appendChild(o);

    o.querySelector('#cwClose').addEventListener('click', close);
    const size = o.querySelector('#cwSize'), sizeVal = o.querySelector('#cwSizeVal');
    size.addEventListener('input', () => { sizeVal.textContent = size.value + ' nm'; });
    size.addEventListener('change', () => {
      const note = o.querySelector('#cwCountNote');
      if (note) note.innerHTML = '<b>配置を計算しています…</b>';
      running = false;
      // 1フレーム描かせてから構築する（大きい箱では数百ms かかるため）
      requestAnimationFrame(() => requestAnimationFrame(() => {
        rebuild(Number(size.value));
        running = true; lastFrame = 0;
        raf = requestAnimationFrame(frame);
      }));
    });
    const speed = o.querySelector('#cwSpeed'), speedVal = o.querySelector('#cwSpeedVal');
    speed.addEventListener('input', () => {
      speedSubsteps = Number(speed.value);
      speedVal.textContent = '×' + speed.value;
    });
    o.querySelector('#cwCrowd').addEventListener('change', e => {
      sim.crowding = e.target.checked;
      sim.resetMSD();
    });
    o.querySelector('#cwResetMsd').addEventListener('click', () => sim.resetMSD());
    document.addEventListener('keydown', e => {
      if (e.key === 'Escape' && overlay && overlay.style.display === 'flex') close();
    });
    return o;
  }

  // 細胞のどこを切り出しているかを示す図
  function drawInset(){
    const host = overlay.querySelector('#cwInset');
    const L = sim ? sim.L : 150;
    // 大腸菌 ≈ 長さ2µm × 直径1µm。図は 2000nm × 1000nm を 180×90 で描く。
    const px = 180 / 2000, boxPx = Math.max(3, L * px);
    host.innerHTML =
      '<svg viewBox="0 0 190 100">' +
        '<rect x="5" y="20" width="180" height="60" rx="30" fill="#12324a" stroke="#2f6d92"/>' +
        '<path d="M30 40 q20 -10 40 0 t40 0 t40 0" fill="none" stroke="#3f8fb5" stroke-width="2" opacity=".7"/>' +
        '<path d="M30 62 q20 10 40 0 t40 0 t40 0" fill="none" stroke="#3f8fb5" stroke-width="2" opacity=".7"/>' +
        '<rect x="' + (95 - boxPx/2) + '" y="' + (50 - boxPx/2) + '" width="' + boxPx + '" height="' + boxPx +
          '" fill="#ffd166" fill-opacity=".85" stroke="#fff3c4"/>' +
        '<text x="95" y="93" text-anchor="middle" fill="#8fb9d0" font-size="9">大腸菌 ≈ 2 µm（黄色が計算範囲）</text>' +
      '</svg>';
    const note = overlay.querySelector('#cwCountNote');
    if (note && sim) {
      note.innerHTML = '一辺 ' + L + ' nm の立方体 = 細胞体積の約 ' +
        (Math.pow(L/1000, 3) / 1.0 * 100).toFixed(2) + ' %。分子数 <b>' + sim.N.toLocaleString() + '</b> 個。';
    }
  }

  window.Crowding = { open: open, close: close };
})();
