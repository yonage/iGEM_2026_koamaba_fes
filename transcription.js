/*
 * 転写・翻訳の3D分子アニメーション
 *
 * 組んだ遺伝子回路をそのままDNAの二重らせんとして並べ、
 *   ① RNAポリメラーゼがプロモーターを探して結合する
 *   ② 二重らせんがほどけて「転写バブル」ができる
 *   ③ 周囲を飛び回るヌクレオチド(NTP)が1つずつ取り込まれ、mRNAが伸びる
 *   ④ ターミネーターでRNAがヘアピンを作り、ポリメラーゼが外れる
 *   ⑤ リボソームがmRNAを3塩基ずつ読み、アミノ酸の鎖(タンパク質)を作る
 *   ⑥ タンパク質の鎖が折りたたまれて放出される
 * の流れを、分子1つ1つが熱でゆれ動く様子とともに描きます。
 *
 * 物理として計算している部分:
 *   - NTP・リン酸・遊離したRNA鎖・タンパク質鎖のブラウン運動
 *     （過減衰ランジュバン方程式: 変位 = 力×移動度×dt + √(2Ddt)·ξ）
 *   - RNA鎖・タンパク質鎖はビーズとばねの高分子モデル（排除体積つき）
 *   - ターミネーター配列の相補塩基どうしの引力 → ヘアピンが自然に折れる
 *   - タンパク質鎖のドメイン内引力 → 球状に折りたたまれる
 * 演出として決めている部分:
 *   - ポリメラーゼやリボソームの移動の向き・速さ、各段階の切り替え
 *
 * 結合のしやすさ（RNAポリメラーゼがプロモーターに結合できる確率）は、
 * ODEシミュレーションで計算した転写速度 k_tx から決めています。
 * 長さの単位は Å（オングストローム）です。
 */
(function(){
  'use strict';
  if (typeof THREE === 'undefined' || !window.MolLook) return;
  const ML = window.MolLook;
  const V3 = THREE.Vector3;

  // ---- DNA の形（B型DNA）----
  const RISE = 3.4, TWIST = 2 * Math.PI / 10.5, GROOVE = 2.35;
  const R_P = 9.0, R_S = 6.3, R_BASE = 2.75;
  // ブロック1つを何塩基対で表すか
  const BP = { promoter: 22, gene: 54, linker: 9, seq: 10, terminator: 16, pad5: 8, pad3: 14 };

  const COL = {
    pad: 0x7b8794, sugar: 0xe9e3d2, rna: 0xf5a623, phos: 0xffb04a,
    rnapA: 0x34609e, rnapB: 0x4f82c4, rnapC: 0x3f948c, sigma: 0xe3c565,
    ssu: 0xb98fd6, lsu: 0x7d4fae, repressor: 0xff7043, activator: 0x26c6da,
    signal: 0xffe066, linkerAA: 0x94a3b8, trna: 0x6fcfb0
  };
  const BASE_COL = { A: 0xff6b6b, T: 0xffd43b, U: 0xf783ac, G: 0x69db7c, C: 0x4dabf7 };
  const SHAPE_COL = { promoter: 0x4c97ff, visible: 0x9966ff, control: 0xf9a825, meta: 0xff8a65, terminator: 0xff6680, seqrecomb: 0x69f0ae };
  const COMP = { A: 'T', T: 'A', G: 'C', C: 'G' };
  const N_NTP = 48, MAX_PPI = 14, MAX_AA = 420, MAX_SIG = 16;

  function proteinColor(part){
    if (part.id === 'enzFast') return 0xf97316;
    if (part.id === 'enzStable') return 0x14b8a6;
    if (part.id === 'DegradingEnzyme') return 0xf59e0b;
    return ({ 'gene-visible': 0x39d98a, 'gene-repressor': 0xff7043, 'gene-activator': 0x26c6da,
              'gene-kill': 0xe53935, 'gene-recomb': 0xa1887f })[part.kind] || 0xf97316;
  }

  // ---- 小道具 ----
  function rng(seed){
    let a = (seed >>> 0) || 1;
    return function(){
      a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function gauss(){
    let u = 0, v = 0;
    while (u === 0) u = Math.random();
    while (v === 0) v = Math.random();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }
  function smooth(e0, e1, x){ const t = Math.max(0, Math.min(1, (x - e0) / (e1 - e0))); return t * t * (3 - 2 * t); }
  function randUnit(out){ out.set(gauss(), gauss(), gauss()); return out.normalize(); }
  function jiggle(v, a){ v.x += gauss() * a; v.y += gauss() * a; v.z += gauss() * a; return v; }

  const colCache = new Map();
  function linCol(hex){ let c = colCache.get(hex); if (!c) { c = ML.lin(hex); colCache.set(hex, c); } return c; }

  function makeMesh(scene, count, detail, mat){
    const mesh = new THREE.InstancedMesh(new THREE.IcosahedronGeometry(1, detail), mat, count);
    mesh.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(count * 3), 3);
    mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    mesh.frustumCulled = false;
    scene.add(mesh);
    return mesh;
  }
  function setInst(mesh, i, x, y, z, s){
    const a = mesh.instanceMatrix.array, o = i * 16;
    a[o] = s; a[o+5] = s; a[o+10] = s; a[o+15] = 1;
    a[o+12] = x; a[o+13] = y; a[o+14] = z;
  }
  function setCol(mesh, i, hex){ const c = linCol(hex); mesh.instanceColor.setXYZ(i, c.r, c.g, c.b); }
  function hideFrom(mesh, from){ for (let i = from; i < mesh.count; i++) setInst(mesh, i, 0, -1e5, 0, 0); }
  function flush(mesh){ mesh.instanceMatrix.needsUpdate = true; if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true; }

  // =====================================================================
  // 回路 → DNA配列
  // =====================================================================
  function buildTrack(blocks, seed){
    const R = rng(seed);
    const regions = [];
    let n = 0;
    const add = (type, part, len) => {
      regions.push({ type, part, start: n, end: n + len,
        color: type === 'pad' ? COL.pad : (SHAPE_COL[part.shape] || COL.pad),
        label: part ? part.label : '', short: part ? (part.short || part.label) : '' });
      n += len;
    };
    add('pad', null, BP.pad5);
    blocks.forEach(part => {
      if (part.cat === 'promoter') add('promoter', part, BP.promoter);
      else if (part.kind === 'terminator') add('terminator', part, BP.terminator);
      else if (part.cat === 'gene') add('gene', part, BP.gene);
      else if (part.kind === 'seq-linker') add('linker', part, BP.linker);
      else add('seq', part, BP.seq);
    });
    add('pad', null, BP.pad3);

    const seq = [];
    for (let i = 0; i < n; i++) seq.push('ATGC'[Math.floor(R() * 4)]);

    const prom = regions.find(r => r.type === 'promoter') || regions[0];
    const tStart = prom.end;
    let tEnd = n, hasTerm = false, termRegion = null;
    for (const r of regions) {
      if (r.start >= tStart && r.type === 'terminator') { tEnd = r.end; hasTerm = true; termRegion = r; break; }
    }

    // 翻訳される読み枠(ORF)。遺伝子-リンカー-遺伝子 は1本の融合タンパク質になる
    const orfs = [];
    const tr = regions.filter(r => r.start >= tStart && r.end <= tEnd);
    let cur = null;
    for (let i = 0; i < tr.length; i++) {
      const r = tr[i], prev = tr[i - 1];
      if (r.type === 'gene') {
        if (cur && prev && prev.type === 'linker' && cur.segs[cur.segs.length - 1] === prev) cur.segs.push(r);
        else { cur = { segs: [r] }; orfs.push(cur); }
      } else if (r.type === 'linker' && cur && prev === cur.segs[cur.segs.length - 1] && tr[i + 1] && tr[i + 1].type === 'gene') {
        cur.segs.push(r);
      } else cur = null;
    }
    orfs.forEach(o => {
      o.start = o.segs[0].start; o.end = o.segs[o.segs.length - 1].end;
      const genes = o.segs.filter(s => s.type === 'gene');
      o.fusion = genes.length > 1;
      o.label = genes.map(g => g.short).join('＋') + (o.fusion ? ' 融合タンパク質' : '');
      // 開始コドン ATG と終止コドン TAA
      seq[o.start] = 'A'; seq[o.start + 1] = 'T'; seq[o.start + 2] = 'G';
      seq[o.end - 3] = 'T'; seq[o.end - 2] = 'A'; seq[o.end - 1] = 'A';
      o.aa = [];
      let dom = 0;
      o.segs.forEach((s, si) => {
        if (s.type === 'gene' && si > 0) dom++;
        for (let b = s.start; b + 3 <= s.end; b += 3) {
          if (b === o.end - 3) break;                       // 終止コドンはアミノ酸にならない
          o.aa.push({
            color: s.type === 'gene' ? proteinColor(s.part) : COL.linkerAA,
            glow: s.type === 'gene' && s.part.kind === 'gene-visible',
            dom: s.type === 'gene' ? dom : -1
          });
        }
      });
    });
    return { regions, seq, n, tStart, tEnd, hasTerm, termRegion, orfs, prom };
  }

  // =====================================================================
  // でこぼこした大きな分子（RNAポリメラーゼ・リボソームなど）
  // =====================================================================
  function blobSpheres(R, inside, lo, hi, minSep, maxN, colorOf){
    const pts = [];
    let tries = 0;
    const p = new V3();
    while (pts.length < maxN && tries < maxN * 90) {
      tries++;
      p.set(lo.x + R() * (hi.x - lo.x), lo.y + R() * (hi.y - lo.y), lo.z + R() * (hi.z - lo.z));
      if (!inside(p)) continue;
      let ok = true;
      for (const q of pts) { if (q.p.distanceToSquared(p) < minSep * minSep) { ok = false; break; } }
      if (!ok) continue;
      const shade = 0.9 + R() * 0.2;
      const c = new THREE.Color(colorOf(p)).multiplyScalar(shade);
      pts.push({ p: p.clone(), r: minSep * (0.64 + R() * 0.2), hex: c.getHex() });
    }
    return pts;
  }

  function Body(scene, spheres, mat){
    this.sp = spheres;
    this.pos = new V3(); this.quat = new THREE.Quaternion(); this.scale = 1; this.visible = false;
    this.wob = new V3(); this.wobR = new V3();
    this.jit = spheres.map(() => new V3());
    this.mesh = makeMesh(scene, Math.max(1, spheres.length), 2, mat);
    spheres.forEach((s, i) => setCol(this.mesh, i, s.hex));
    this._q = new THREE.Quaternion(); this._e = new THREE.Euler(); this._v = new V3();
  }
  const OU_TH = 0.06;
  Body.prototype.update = function(h){
    const sh = Math.sqrt(h);
    const ou = (v, amp) => {
      const sg = amp * Math.sqrt(2 * OU_TH);
      v.x += -v.x * OU_TH * h + gauss() * sg * sh;
      v.y += -v.y * OU_TH * h + gauss() * sg * sh;
      v.z += -v.z * OU_TH * h + gauss() * sg * sh;
    };
    ou(this.wob, 1.1); ou(this.wobR, 0.025);
    this.jit.forEach(j => ou(j, 0.45));
  };
  Body.prototype.frameQuat = function(){
    this._e.set(this.wobR.x, this.wobR.y, this.wobR.z);
    return this._q.setFromEuler(this._e).premultiply(this.quat);
  };
  // 局所座標 → 世界座標（ゆれは位置・向きのみ反映、粒子ごとの細かいゆれは除く）
  Body.prototype.toWorld = function(local, out){
    const q = this.frameQuat();
    return out.copy(local).applyQuaternion(q).multiplyScalar(this.scale).add(this.pos).add(this.wob);
  };
  Body.prototype.render = function(){
    if (!this.visible || this.scale < 0.01) { this.mesh.visible = false; return; }
    this.mesh.visible = true;
    const q = this.frameQuat(), v = this._v;
    for (let i = 0; i < this.sp.length; i++) {
      const s = this.sp[i];
      v.copy(s.p).add(this.jit[i]).applyQuaternion(q).multiplyScalar(this.scale).add(this.pos).add(this.wob);
      setInst(this.mesh, i, v.x, v.y, v.z, s.r * this.scale);
    }
    flush(this.mesh);
  };
  Body.prototype.dispose = function(scene){ scene.remove(this.mesh); this.mesh.geometry.dispose(); };

  // RNAポリメラーゼ: DNAを抱え込むクランプ形。+z側にすき間が開いていて
  // 中の転写バブルと、NTPが入ってくる様子が見える。
  function rnapSpheres(){
    const R = rng(7);
    return blobSpheres(R, p => {
      const e = (p.x/42)**2 + (p.y/40)**2 + (p.z/38)**2;
      if (e > 1 || e < 0.42) return false;
      if (p.y*p.y + p.z*p.z < 15.5*15.5) return false;              // DNAの通り道
      if (p.z > 13 && p.y > -17 && p.y < 23) return false;           // 正面のすき間
      // RNAの出口トンネル
      const t = Math.max(0, Math.min(1, (p.y - 20) / 28));
      const cx = -26 + 2 * t, cz = 6;
      if (p.y > 16 && Math.hypot(p.x - cx, p.z - cz) < 8.5) return false;
      return true;
    }, new V3(-42,-40,-38), new V3(42,40,38), 8.6, 260,
    p => p.x > 26 ? COL.rnapC : (p.y > 4 ? COL.rnapB : COL.rnapA));
  }
  function sigmaSpheres(){
    const R = rng(11);
    return blobSpheres(R, p => ((p.x+34)/14)**2 + ((p.y+22)/11)**2 + ((p.z-16)/12)**2 <= 1,
      new V3(-48,-33,4), new V3(-20,-11,28), 6.2, 30, () => COL.sigma);
  }
  function ssuSpheres(){
    const R = rng(21);
    return blobSpheres(R, p => {
      const e = (p.x/46)**2 + ((p.y+20)/20)**2 + (p.z/40)**2;
      return e <= 1 && e >= 0.3 && p.y < -5;
    }, new V3(-46,-40,-40), new V3(46,0,40), 9.5, 170, () => COL.ssu);
  }
  function lsuSpheres(){
    const R = rng(23);
    return blobSpheres(R, p => {
      const e = (p.x/54)**2 + ((p.y-32)/34)**2 + (p.z/50)**2;
      if (e > 1 || e < 0.3 || p.y < 6) return false;
      if (Math.hypot(p.x, p.z) < 7 && p.y > 30) return false;        // ペプチドの出口トンネル
      return true;
    }, new V3(-54,6,-50), new V3(54,66,50), 10.5, 240, () => COL.lsu);
  }
  // tRNA: L字形。下の端(アンチコドン)でmRNAを読み、もう一方の端にアミノ酸をのせて運ぶ
  function trnaSpheres(){
    const pts = [];
    for (let i = 0; i <= 4; i++) pts.push({ p: new V3(0, -i * 4.4, 0), r: 3.1, hex: COL.trna });
    for (let i = 1; i <= 3; i++) pts.push({ p: new V3(i * 4.4, 1.5 * i, 0), r: 3.1, hex: COL.trna });
    pts.push({ p: new V3(4 * 4.4 + 1, 6.5, 0), r: 3.8, hex: 0xffffff });   // 運んでいるアミノ酸
    return pts;
  }
  const TRNA_AA = 8;

  function factorSpheres(hex, seed){
    const R = rng(seed);
    return blobSpheres(R, p => (p.x/17)**2 + (p.y/12)**2 + (p.z/13)**2 <= 1,
      new V3(-17,-12,-13), new V3(17,12,13), 6.4, 34, () => hex);
  }

  // =====================================================================
  // ビーズとばねの高分子（RNA鎖・タンパク質鎖）
  // =====================================================================
  function Chain(o){
    this.b = [];
    this.b0 = o.b0; this.k2r = o.k2r; this.k2 = o.k2; this.excl = o.excl;
    this.collapse = o.collapse || 0; this.noise = o.noise; this.hasBase = !!o.hasBase;
    this.pairs = []; this.scale = 1; this.drift = null; this.age = 0; this.dying = false;
  }
  Chain.prototype.add = function(p, extra){
    const bd = Object.assign({ p: p.clone(), q: p.clone(), u: randUnit(new V3()), f: new V3(),
      pin: new V3(), pinned: false, pinQ: new V3(), pinnedQ: false }, extra);
    this.b.push(bd);
    return bd;
  };
  Chain.prototype.unpinAll = function(){ this.b.forEach(x => { x.pinned = false; x.pinnedQ = false; }); };
  const _t = new V3();
  Chain.prototype.step = function(h, env){
    const b = this.b, n = b.length;
    if (!n) return;
    for (let i = 0; i < n; i++) b[i].f.set(0, 0, 0);
    const spring = (i, j, r0, k) => {
      _t.subVectors(b[j].p, b[i].p);
      const r = _t.length() || 1e-6, m = k * (r - r0) / r;
      b[i].f.addScaledVector(_t, m); b[j].f.addScaledVector(_t, -m);
    };
    for (let i = 0; i < n - 1; i++) spring(i, i + 1, this.b0, 1);
    for (let i = 0; i < n - 2; i++) spring(i, i + 2, this.k2r, this.k2);
    for (const pr of this.pairs) {
      if (pr[1] < n && !b[pr[0]].pinned && !b[pr[1]].pinned) spring(pr[0], pr[1], 6.4, 0.4);
    }
    const ex = this.excl, ex2 = ex * ex, cr2 = (ex * 2.6) ** 2;
    for (let i = 0; i < n; i++) {
      const pi = b[i].p;
      for (let j = i + 3; j < n; j++) {
        const pj = b[j].p;
        const dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
        const d2 = dx*dx + dy*dy + dz*dz;
        let m = 0;
        if (d2 < ex2) { const d = Math.sqrt(d2) || 1e-6; m = -0.6 * (ex - d) / d; }
        else if (this.collapse && d2 < cr2 && b[i].dom >= 0 && b[i].dom === b[j].dom) {
          const d = Math.sqrt(d2); m = this.collapse * (d - this.excl * 1.15) / d;
        }
        if (m) {
          b[i].f.x += dx*m; b[i].f.y += dy*m; b[i].f.z += dz*m;
          b[j].f.x -= dx*m; b[j].f.y -= dy*m; b[j].f.z -= dz*m;
        }
      }
    }
    const sh = Math.sqrt(h), nz = this.noise;
    for (let i = 0; i < n; i++) {
      const x = b[i];
      if (x.pinned) { x.p.lerp(x.pin, 1 - Math.exp(-0.3 * h)); }
      else {
        if (env) {
          for (const s of env.spheres) {
            _t.subVectors(x.p, s.c);
            const d = _t.length();
            if (d < s.r) x.f.addScaledVector(_t, (s.r - d) / (d || 1e-6) * 0.7);
          }
          if (x.p.x > env.x0 && x.p.x < env.x1) {
            const ryz = Math.hypot(x.p.y, x.p.z) || 1e-6;
            if (ryz < 11.5) { const m = (11.5 - ryz) / ryz * 0.7; x.f.y += x.p.y * m; x.f.z += x.p.z * m; }
          }
        }
        if (this.drift) x.f.add(this.drift);
        x.f.clampLength(0, 30);
        x.p.addScaledVector(x.f, 0.12 * h);
        x.p.x += gauss() * nz * sh; x.p.y += gauss() * nz * sh; x.p.z += gauss() * nz * sh;
      }
      if (this.hasBase) {
        if (x.pinnedQ) x.q.lerp(x.pinQ, 1 - Math.exp(-0.3 * h));
        else {
          jiggle(x.u, 0.12 * sh).normalize();
          _t.copy(x.p).addScaledVector(x.u, 4.2);
          x.q.lerp(_t, 1 - Math.exp(-0.4 * h));
        }
      }
    }
    this.age += h / 60;
  };

  // =====================================================================
  // ビュー
  // =====================================================================
  let renderer = null, scene = null, camera = null, backdrop = null;
  let host = null, mountEl = null, ro = null, raf = null, running = false, lastT = 0;
  let W = null;                         // いま表示している「世界」
  let speed = 1, autoCam = true, userOrbit = { theta: Math.PI / 2 - 0.45, phi: 1.2, dist: 220 };
  let dragging = false, lastPtr = null, camT = 0;
  const focus = new V3(), camGoal = new V3();
  let ui = {};

  function initRenderer(){
    renderer = new THREE.WebGLRenderer({ antialias: true });
    ML.setupRenderer(renderer);
    scene = new THREE.Scene();
    ML.setupScene(scene, 0.0026);
    camera = new THREE.PerspectiveCamera(40, 1, 1, 3000);
    backdrop = new ML.Backdrop(scene, 46, 420);
    const dom = renderer.domElement;
    dom.addEventListener('pointerdown', e => { dragging = true; lastPtr = { x: e.clientX, y: e.clientY }; });
    window.addEventListener('pointermove', e => {
      if (!dragging) return;
      if (autoCam) { autoCam = false; syncButtons(); }
      userOrbit.theta -= (e.clientX - lastPtr.x) * 0.007;
      userOrbit.phi = Math.max(0.15, Math.min(Math.PI - 0.15, userOrbit.phi - (e.clientY - lastPtr.y) * 0.007));
      lastPtr = { x: e.clientX, y: e.clientY };
    });
    window.addEventListener('pointerup', () => { dragging = false; });
    dom.addEventListener('wheel', e => {
      e.preventDefault();
      userOrbit.dist = Math.max(70, Math.min(700, userOrbit.dist * (1 + e.deltaY * 0.0012)));
    }, { passive: false });
  }

  // ---------------------------------------------------------------------
  function buildWorld(circuit){
    disposeWorld();
    const tr = buildTrack(circuit.blocks, hashStr(circuit.blocks.map(b => b.id).join('|')));
    const w = {
      circuit, tr, n: tr.n,
      pBind: Math.max(0, Math.min(1, (circuit.k_tx || 0) / 1.5)),
      phase: 'search', t: 0, s: tr.tStart - 1, sF: tr.tStart - 1, bubble: 0,
      rna: null, prot: null, products: [], ntps: [], ppis: [], sigs: [],
      cap: null, capT: 0, trnaIn: null, mRNAs: 0, proteins: 0, binds: 0, fails: 0,
      kF: 0, kGoal: 0, orfIdx: 0, aaIdx: 0, codonT: 0,
      sigmaOff: false, sigmaVel: new V3(), factorAway: 0
    };
    w.xOf = i => (i - tr.tStart) * RISE;
    W = w;

    // DNA
    const mDNA = ML.material(0xffffff, { vertexColor: true, roughness: 0.55, rimStr: 0.5 });
    w.mats = [mDNA];
    w.dP = makeMesh(scene, w.n * 2, 2, mDNA);
    w.dS = makeMesh(scene, w.n * 2, 1, mDNA);
    w.dB = makeMesh(scene, w.n * 2, 1, mDNA);
    w.helix = new Float32Array(w.n * 18);
    for (let i = 0; i < w.n; i++) {
      const reg = tr.regions.find(r => i >= r.start && i < r.end);
      const c = new THREE.Color(reg.color);
      setCol(w.dP, 2*i, c.getHex());
      setCol(w.dP, 2*i+1, c.clone().multiplyScalar(0.72).getHex());
      setCol(w.dS, 2*i, c.clone().lerp(new THREE.Color(0xffffff), 0.45).getHex());
      setCol(w.dS, 2*i+1, c.clone().lerp(new THREE.Color(0xffffff), 0.3).getHex());
      setCol(w.dB, 2*i, BASE_COL[tr.seq[i]]);
      setCol(w.dB, 2*i+1, BASE_COL[COMP[tr.seq[i]]]);
      const th = i * TWIST, ph = th + GROOVE, x = w.xOf(i), h = w.helix, o = i * 18;
      // 鎖A(コード鎖): リン酸・糖・塩基
      h[o] = x;       h[o+1] = R_P*Math.cos(th);        h[o+2] = R_P*Math.sin(th);
      h[o+3] = x+0.8; h[o+4] = R_S*Math.cos(th+0.32);   h[o+5] = R_S*Math.sin(th+0.32);
      h[o+6] = x+0.8; h[o+7] = R_BASE*Math.cos(th+0.32); h[o+8] = R_BASE*Math.sin(th+0.32);
      // 鎖B(鋳型鎖)
      h[o+9] = x;      h[o+10] = R_P*Math.cos(ph);        h[o+11] = R_P*Math.sin(ph);
      h[o+12] = x-0.8; h[o+13] = R_S*Math.cos(ph-0.32);   h[o+14] = R_S*Math.sin(ph-0.32);
      h[o+15] = x-0.8; h[o+16] = R_BASE*Math.cos(ph-0.32); h[o+17] = R_BASE*Math.sin(ph-0.32);
    }

    // RNA
    const nRNA = Math.max(1, tr.tEnd - tr.tStart);
    const mRNA = ML.material(0xffffff, { vertexColor: true, roughness: 0.5, rimStr: 0.55 });
    w.mats.push(mRNA);
    w.rP = makeMesh(scene, nRNA, 1, mRNA);
    w.rB = makeMesh(scene, nRNA, 1, mRNA);
    for (let j = 0; j < nRNA; j++) {
      setCol(w.rP, j, COL.rna);
      const b = tr.seq[tr.tStart + j];
      setCol(w.rB, j, BASE_COL[b === 'T' ? 'U' : b]);
    }

    // NTP（塩基・糖・三リン酸）とピロリン酸
    const mNTP = ML.material(0xffffff, { vertexColor: true, roughness: 0.45, rimStr: 0.6 });
    w.mats.push(mNTP);
    w.nB = makeMesh(scene, N_NTP, 1, mNTP);
    w.nS = makeMesh(scene, N_NTP, 1, mNTP);
    w.nP = makeMesh(scene, N_NTP * 3 + MAX_PPI * 2, 1, mNTP);
    for (let i = 0; i < N_NTP * 3 + MAX_PPI * 2; i++) setCol(w.nP, i, COL.phos);
    for (let i = 0; i < N_NTP; i++) {
      setCol(w.nS, i, COL.sugar);
      const nt = { c: new V3(), q: new THREE.Quaternion(), type: 'AUGC'[i % 4], state: 'free' };
      w.ntps.push(nt);
      respawnNTP(nt, new V3(0, 0, 0), true);
    }

    // シグナル分子（シグナル応答プロモーターのときだけ）
    const sig = circuit.signal || 0;
    w.sigMesh = makeMesh(scene, MAX_SIG * 2, 1, mNTP);
    if (/^prom-signal/.test(circuit.promKind) && sig > 0.02) {
      const nSig = Math.min(MAX_SIG, 2 + Math.round(sig * 14));
      const pc = new V3(w.xOf((tr.prom.start + tr.prom.end) / 2), 0, 0);
      for (let i = 0; i < nSig; i++) {
        w.sigs.push({ c: pc.clone().add(randUnit(new V3()).multiplyScalar(20 + Math.random() * 50)), q: new THREE.Quaternion(), home: pc });
      }
    }
    for (let i = 0; i < MAX_SIG * 2; i++) setCol(w.sigMesh, i, COL.signal);

    // アミノ酸（GFPは光らせる）
    const mAA = ML.material(0xffffff, { vertexColor: true, roughness: 0.5, rimStr: 0.55 });
    const mGlow = ML.material(0xffffff, { vertexColor: true, roughness: 0.4, emissive: 0x16b85a, emissiveIntensity: 1.5, rim: 0xa8ffd0, rimStr: 0.85 });
    w.mats.push(mAA, mGlow);
    w.aa = makeMesh(scene, MAX_AA, 1, mAA);
    w.aaG = makeMesh(scene, MAX_AA, 1, mGlow);

    // 大きな分子
    const mBody = ML.material(0xffffff, { vertexColor: true, roughness: 0.72, rimStr: 0.45, rimPow: 2.8 });
    w.mats.push(mBody);
    w.rnap = new Body(scene, rnapSpheres(), mBody);
    w.sigma = new Body(scene, sigmaSpheres(), mBody);
    w.ssu = new Body(scene, ssuSpheres(), mBody);
    w.lsu = new Body(scene, lsuSpheres(), mBody);
    w.trnas = [0, 1, 2].map(() => { const t = new Body(scene, trnaSpheres(), mBody); t.state = 'idle'; t.t = 0; return t; });
    w.trnaNext = 0;
    w.factor = null;
    if (circuit.promKind === 'prom-repressor' && w.pBind < 0.6) {
      w.factor = new Body(scene, factorSpheres(COL.repressor, 31), mBody);
      w.factor.kind = 'repressor';
      w.factor.home = new V3(w.xOf((tr.prom.start + tr.prom.end) / 2), 6, 16);   // プロモーターの手前側をふさぐ
    } else if (circuit.promKind === 'prom-activator' && w.pBind > 0.15) {
      w.factor = new Body(scene, factorSpheres(COL.activator, 33), mBody);
      w.factor.kind = 'activator';
      w.factor.home = new V3(w.xOf(tr.prom.start + 3), -15, 4);
    }
    if (w.factor) { w.factor.visible = true; w.factor.pos.copy(w.factor.home); }

    w.env = { spheres: [], x0: w.xOf(0) - 2, x1: w.xOf(w.n - 1) + 2 };
    buildMap();
    startPhase('search');
    focusTarget(focus);
  }

  function disposeWorld(){
    if (!W) return;
    ['dP','dS','dB','rP','rB','nB','nS','nP','aa','aaG','sigMesh'].forEach(k => {
      if (W[k]) { scene.remove(W[k]); W[k].geometry.dispose(); }
    });
    ['rnap','sigma','ssu','lsu','factor'].forEach(k => { if (W[k]) W[k].dispose(scene); });
    W.trnas.forEach(t => t.dispose(scene));
    W.mats.forEach(m => m.dispose());
    W = null;
  }

  function hashStr(s){
    let h = 2166136261;
    for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 16777619); }
    return h >>> 0;
  }

  // ---------------------------------------------------------------------
  function respawnNTP(nt, center, anywhere){
    const v = randUnit(new V3()).multiplyScalar(anywhere ? 30 + Math.random() * 80 : 85 + Math.random() * 20);
    nt.c.copy(center).add(v);
    nt.state = 'free';
    randQuat(nt.q);
  }
  function randQuat(q){ q.set(gauss(), gauss(), gauss(), gauss()).normalize(); return q; }
  const _dq = new THREE.Quaternion(), _ax = new V3();
  function tumble(q, amt){
    randUnit(_ax); _dq.setFromAxisAngle(_ax, gauss() * amt); q.multiply(_dq).normalize();
  }

  // 局所座標の点を各分子の世界座標へ
  const L = (x, y, z) => new V3(x, y, z);
  const LOC = {
    window: L(10, 3, 34), exitCtl: L(-26, 22, 8), exit: L(-24, 46, 6),
    tunnel: L(0, 72, 0), sigmaAway: L(-60, -40, 50)
  };
  const tmpA = new V3(), tmpB = new V3(), tmpC = new V3(), tmpD = new V3();

  // 転写バブルの開き具合（0=二重らせん, 1=完全に開いた状態）
  function openAt(i){
    if (W.bubble <= 0) return 0;
    const d = i - W.sF;
    return W.bubble * smooth(-12, -9, d) * (1 - smooth(2, 4.5, d));
  }

  // RNA:DNAハイブリッド内のRNA塩基・骨格の位置（DNAと同じ座標系）
  function hybridBase(b, out){ return out.set(W.xOf(b), -2.6, 3.6); }
  function hybridBack(b, out){ return out.set(W.xOf(b), 1.8, 6.8); }

  // =====================================================================
  // 段階（フェーズ）の進行
  // =====================================================================
  function startPhase(ph){
    const w = W;
    w.phase = ph; w.t = 0;
    if (ph === 'search') {
      w.success = Math.random() < w.pBind;
      w.rnap.visible = true; w.rnap.scale = 1;
      w.rnap.pos.set((Math.random() - 0.5) * 80, 30 + Math.random() * 40, 170);
      w.rnap.quat.identity();
      w.sigma.visible = true; w.sigma.scale = 1; w.sigmaOff = false;
      w.s = w.tr.tStart - 1; w.sF = w.s; w.bubble = 0;
      // 結合できないときは、リプレッサーなどにぶつかる手前で止まって離れていく
      w.goal = w.success ? new V3(w.xOf(w.s), 0, 0) : new V3(w.xOf(w.s) - 10, 74, 36);
    } else if (ph === 'bind') {
      w.binds++;
    } else if (ph === 'bounce') {
      w.fails++;
      w.goal = w.rnap.pos.clone().add(new V3((Math.random() - 0.5) * 120, 50, 150));
    } else if (ph === 'elong') {
      w.rna = new Chain({ b0: 6.0, k2r: 10.5, k2: 0.15, excl: 5.0, noise: 0.75, hasBase: true });
      if (w.tr.hasTerm) {
        const a = w.tr.termRegion.start - w.tr.tStart;
        for (let m = 0; m < 5; m++) w.rna.pairs.push([a + 1 + m, a + 14 - m]);
      }
      w.cap = null;
    } else if (ph === 'release') {
      w.mRNAs++;
      w.goal = w.rnap.pos.clone().add(new V3(0, 60, 180));
      if (w.rna) w.rna.unpinAll();
    } else if (ph === 'riboIn') {
      w.orfIdx = 0;
      const j0 = w.tr.orfs[0].start - w.tr.tStart;
      const b = w.rna.b[Math.min(j0, w.rna.b.length - 1)];
      w.ssu.visible = true; w.ssu.scale = 1; w.ssu.quat.identity();
      w.ssu.pos.copy(b.p).add(new V3((Math.random() - 0.5) * 60, -70, 130));
      w.lsu.visible = false;
      w.kF = j0; w.kGoal = j0;
    } else if (ph === 'lsuIn') {
      w.lsu.visible = true; w.lsu.scale = 1; w.lsu.quat.identity();
      w.lsu.pos.copy(w.ssu.pos).add(new V3(20, 150, 60));
    } else if (ph === 'translate') {
      newProtein();
      w.codonT = 0;
    } else if (ph === 'riboOut') {
      w.rna.unpinAll();
      if (w.trnaIn) { w.trnaIn.state = 'out'; w.trnaIn.t = 0; w.trnaIn = null; }
    } else if (ph === 'fade') {
      w.rna && (w.rna.dying = true);
    }
  }

  function newProtein(){
    const w = W, orf = w.tr.orfs[w.orfIdx];
    w.prot = new Chain({ b0: 5.4, k2r: 8.8, k2: 0.2, excl: 6.6, collapse: 0.05, noise: 0.5 });
    w.prot.orf = orf;
    w.aaIdx = 0;
    w.kGoal = orf.start - w.tr.tStart;
  }

  function releaseProtein(){
    const w = W;
    if (!w.prot) return;
    w.prot.unpinAll();
    w.prot.drift = randUnit(new V3()).add(new V3(0, 1.2, 0.4)).normalize().multiplyScalar(1.6);
    w.prot.age = 0;
    w.products.push(w.prot);
    w.proteins++;
    w.prot = null;
  }

  function nearestNTP(type, from){
    let best = null, bd = 1e9;
    for (const nt of W.ntps) {
      if (nt.state !== 'free' || nt.type !== type) continue;
      const d = nt.c.distanceToSquared(from);
      if (d < bd) { bd = d; best = nt; }
    }
    return best;
  }

  function step(dt){
    const w = W, tr = w.tr, h = dt * 60;
    w.t += dt;

    // 大きな分子のゆれ
    ['rnap','sigma','ssu','lsu','factor'].forEach(k => { if (w[k] && w[k].visible) w[k].update(h); });

    // ---- 段階ごとの動き ----
    const ph = w.phase;
    if (ph === 'search') {
      const k = 1 - Math.exp(-0.85 * dt);
      w.rnap.pos.lerp(w.goal, k);
      jiggle(w.rnap.pos, 0.35 * Math.sqrt(h));
      if (w.rnap.pos.distanceTo(w.goal) < 2.2 || w.t > 7) startPhase(w.success ? 'bind' : 'bounce');
    } else if (ph === 'bounce') {
      w.rnap.pos.lerp(w.goal, 1 - Math.exp(-0.7 * dt));
      jiggle(w.rnap.pos, 0.4 * Math.sqrt(h));
      if (w.t > 1.8) w.rnap.scale = Math.max(0, 1 - (w.t - 1.8) / 0.8);
      w.sigma.scale = w.rnap.scale;
      if (w.t > 2.9) startPhase('search');
    } else if (ph === 'bind') {
      w.rnap.pos.lerp(tmpA.set(w.xOf(w.sF), 0, 0), 1 - Math.exp(-3 * dt));
      w.bubble = smooth(0.2, 1.8, w.t);
      if (w.t > 2.1) startPhase('elong');
    } else if (ph === 'elong') {
      stepElongation(dt, h);
    } else if (ph === 'term') {
      w.rnap.pos.set(w.xOf(w.sF), 0, 0);
      if (w.t > 1.2) startPhase('release');
    } else if (ph === 'release') {
      w.bubble = Math.max(0, 1 - w.t / 1.4);
      w.rnap.pos.lerp(w.goal, 1 - Math.exp(-0.55 * dt));
      jiggle(w.rnap.pos, 0.35 * Math.sqrt(h));
      if (w.t > 2.0) w.rnap.scale = Math.max(0, 1 - (w.t - 2.0) / 0.9);
      if (w.t > 3.0) { w.rnap.visible = false; startPhase(tr.orfs.length ? 'riboIn' : 'fade'); }
    } else if (ph === 'riboIn') {
      const j0 = tr.orfs[0].start - tr.tStart;
      const b = w.rna.b[Math.min(j0, w.rna.b.length - 1)];
      w.ssu.pos.lerp(b.p, 1 - Math.exp(-1.1 * dt));
      jiggle(w.ssu.pos, 0.35 * Math.sqrt(h));
      if (w.ssu.pos.distanceTo(b.p) < 4 || w.t > 4) startPhase('lsuIn');
    } else if (ph === 'lsuIn') {
      w.lsu.pos.lerp(w.ssu.pos, 1 - Math.exp(-2.4 * dt));
      if (w.t > 1.6) { w.lsu.pos.copy(w.ssu.pos); startPhase('translate'); }
    } else if (ph === 'translate') {
      stepTranslation(dt);
    } else if (ph === 'scan') {
      w.kF += Math.sign(w.kGoal - w.kF) * Math.min(Math.abs(w.kGoal - w.kF), dt * 14);
      if (Math.abs(w.kGoal - w.kF) < 0.01) startPhase('translate');
    } else if (ph === 'riboOut') {
      w.ssu.pos.y -= 30 * dt; w.ssu.pos.z += 40 * dt;
      w.lsu.pos.y += 40 * dt; w.lsu.pos.z += 30 * dt;
      const sc = Math.max(0, 1 - Math.max(0, w.t - 1.2) / 1.0);
      w.ssu.scale = w.lsu.scale = sc;
      if (w.t > 2.3) { w.ssu.visible = w.lsu.visible = false; startPhase('fade'); }
    } else if (ph === 'fade') {
      if (w.rna) w.rna.scale = Math.max(0, 1 - w.t / 2.5);
      if (w.t > 2.7) { w.rna = null; startPhase('search'); }
    }

    // リボソームはSSUを基準に、LSUはSSUにくっつく
    if ((ph === 'translate' || ph === 'scan') && w.ssu.visible) {
      jiggle(w.ssu.pos, 0.12 * Math.sqrt(h));
      w.lsu.pos.copy(w.ssu.pos);
    }

    // σ因子: 転写が始まるとポリメラーゼから離れていく
    if (w.sigma.visible) {
      if (!w.sigmaOff) { w.sigma.pos.copy(w.rnap.pos); w.sigma.quat.copy(w.rnap.quat); w.sigma.wob.copy(w.rnap.wob); w.sigma.wobR.copy(w.rnap.wobR); }
      else {
        w.sigma.pos.addScaledVector(w.sigmaVel, dt);
        jiggle(w.sigma.pos, 0.4 * Math.sqrt(h));
        w.sigma.scale = Math.max(0, w.sigma.scale - dt / 2.5);
        if (w.sigma.scale <= 0) w.sigma.visible = false;
      }
    }

    // リプレッサー/アクチベーター: ポリメラーゼが結合する間だけ場所を譲る
    if (w.factor) {
      const yield_ = w.factor.kind === 'repressor' &&
        ((ph === 'search' && w.success) || ph === 'bind' || (ph === 'elong' && w.s - tr.tStart < 14));
      w.factorAway += ((yield_ ? 1 : 0) - w.factorAway) * (1 - Math.exp(-2 * dt));
      w.factor.pos.copy(w.factor.home).add(tmpA.set(-30 * w.factorAway, 80 * w.factorAway, 30 * w.factorAway));
    }

    // 障害物（大きな分子）
    w.env.spheres.length = 0;
    if (w.rnap.visible && w.rnap.scale > 0.5) w.env.spheres.push({ c: w.rnap.pos.clone(), r: 40 });
    if (w.ssu.visible && w.ssu.scale > 0.5) {
      w.env.spheres.push({ c: w.ssu.pos.clone().add(tmpA.set(0, -20, 0)), r: 32 });
      if (w.lsu.visible) w.env.spheres.push({ c: w.lsu.pos.clone().add(tmpA.set(0, 32, 0)), r: 44 });
    }

    // RNA鎖のピン留め（ハイブリッド・出口トンネル・リボソーム内）
    if (w.rna) pinRNA();
    if (w.rna) subStep(w.rna, h);
    if (w.prot) {
      if (w.prot.b.length) {
        const last = w.prot.b[w.prot.b.length - 1];
        last.pinned = true; w.ssu.toWorld(LOC.tunnel, last.pin);
      }
      subStep(w.prot, h);
    }
    for (let i = w.products.length - 1; i >= 0; i--) {
      const p = w.products[i];
      if (p.age > 2.5) p.drift = null;
      if (w.products.length - i > 4 || p.age > 40) p.dying = true;
      if (p.dying) { p.scale -= dt / 1.5; if (p.scale <= 0) { w.products.splice(i, 1); continue; } }
      subStep(p, h);
    }

    stepNTPs(dt, h);
    stepTRNAs(dt, h);
  }

  function subStep(chain, h){
    const ns = Math.max(1, Math.ceil(h / 1.5));
    for (let k = 0; k < ns; k++) chain.step(h / ns, W.env);
  }

  function pinRNA(){
    const w = W, tr = w.tr, b = w.rna.b, ph = w.phase;
    if (ph === 'elong' || ph === 'term' || ph === 'runoff') {
      const h0 = hybridBack(w.s - 8, tmpB);
      const ctl = w.rnap.toWorld(LOC.exitCtl, tmpC);
      const ex = w.rnap.toWorld(LOC.exit, tmpD);
      for (let j = 0; j < b.length; j++) {
        const bb = tr.tStart + j, d = w.s - bb, x = b[j];
        if (d < 8) {
          x.pinned = true; hybridBack(bb, x.pin);
          x.pinnedQ = true; hybridBase(bb, x.pinQ);
        } else if (d < 14) {
          const t = (d - 7) / 6, u = 1 - t;
          x.pin.set(u*u*h0.x + 2*u*t*ctl.x + t*t*ex.x, u*u*h0.y + 2*u*t*ctl.y + t*t*ex.y, u*u*h0.z + 2*u*t*ctl.z + t*t*ex.z);
          x.pinned = true; x.pinnedQ = false;
        } else { x.pinned = false; x.pinnedQ = false; }
      }
    } else if (ph === 'lsuIn' || ph === 'translate' || ph === 'scan') {
      for (let j = 0; j < b.length; j++) {
        const x = b[j], d = j - w.kF;
        if (Math.abs(d) <= 8) {
          x.pinned = true; w.ssu.toWorld(tmpA.set(d * 5.6, 0, 0), x.pin);
          x.pinnedQ = true; w.ssu.toWorld(tmpA.set(d * 5.6, -4.4, 1.0), x.pinQ);
        } else { x.pinned = false; x.pinnedQ = false; }
      }
    }
  }

  function stepElongation(dt, h){
    const w = W, tr = w.tr;
    w.sF += (w.s - w.sF) * (1 - Math.exp(-14 * dt));
    w.rnap.pos.set(w.xOf(w.sF), 0, 0);
    w.bubble = 1;

    const next = w.s + 1;
    const end = tr.hasTerm ? tr.tEnd : w.n;
    if (next >= end) { startPhase(tr.hasTerm ? 'term' : 'release'); return; }

    if (!w.cap) {
      const base = tr.seq[next], need = base === 'T' ? 'U' : base;
      let nt = nearestNTP(need, w.rnap.pos);
      if (!nt || nt.c.distanceTo(w.rnap.pos) > 130) {
        // 近くに合うNTPがなければ、正面のすき間の外から1つ流れてくることにする
        nt = W.ntps.find(x => x.state === 'free') || W.ntps[0];
        nt.c.copy(w.rnap.toWorld(LOC.window, tmpA)).add(randUnit(tmpB).multiplyScalar(25)).add(tmpC.set(0, 0, 35));
        nt.type = need;
      }
      nt.state = 'cap1';
      w.cap = nt; w.capT = 0;
    }
    const nt = w.cap;
    w.capT += dt;
    if (nt.state === 'cap1') {
      const W1 = w.rnap.toWorld(LOC.window, tmpA);
      nt.c.lerp(W1, 1 - Math.exp(-19 * dt));
      jiggle(nt.c, 0.5 * Math.sqrt(h));
      if (nt.c.distanceTo(W1) < 8 || w.capT > 0.6) nt.state = 'cap2';
    } else {
      const F = hybridBase(next, tmpA);
      nt.c.lerp(F, 1 - Math.exp(-21 * dt));
      if (nt.c.distanceTo(F) < 2.2 || w.capT > 1.1) {
        // 取り込み: 骨格に連結し、ピロリン酸(PPi)が外れる
        const bead = w.rna.add(hybridBack(next, tmpB));
        bead.q.copy(F);
        if (w.ppis.length >= MAX_PPI) w.ppis.shift();
        w.ppis.push({ c: F.clone(), q: randQuat(new THREE.Quaternion()), t: 0,
          v: w.rnap.toWorld(LOC.window, tmpC).sub(F).normalize().multiplyScalar(26) });
        respawnNTP(nt, w.rnap.pos, false);
        nt.type = 'AUGC'[Math.floor(Math.random() * 4)];
        w.cap = null;
        w.s = next;
        if (!w.sigmaOff && w.s - tr.tStart >= 5) {
          w.sigmaOff = true;
          w.sigmaVel.copy(LOC.sigmaAway).normalize().multiplyScalar(22);
        }
      }
    }
  }

  function stepTranslation(dt){
    const w = W, orf = w.prot && w.prot.orf;
    if (!orf) return;
    w.kF += (w.kGoal - w.kF) * (1 - Math.exp(-8 * dt));
    w.codonT += dt;
    if (w.aaIdx < orf.aa.length && !w.trnaIn) callTRNA(orf.aa[w.aaIdx]);
    if (w.codonT > 0.62) {
      w.codonT = 0;
      if (w.aaIdx < orf.aa.length) {
        if (w.trnaIn) { w.trnaIn.state = 'out'; w.trnaIn.t = 0; setCol(w.trnaIn.mesh, TRNA_AA, COL.trna); flush(w.trnaIn.mesh); w.trnaIn = null; }
        const a = orf.aa[w.aaIdx++];
        const p0 = w.ssu.toWorld(LOC.tunnel, tmpA);
        if (w.prot.b.length) {
          const prev = w.prot.b[w.prot.b.length - 1];
          prev.pinned = false;
          prev.p.add(tmpB.set(gauss(), 2.5, gauss()));
        }
        w.prot.add(p0, { color: a.color, glow: a.glow, dom: a.dom });
        w.kGoal = (orf.start - w.tr.tStart) + w.aaIdx * 3;
      } else {
        w.kGoal = (orf.end - w.tr.tStart);          // 終止コドンを読む
        releaseProtein();
        w.orfIdx++;
        if (w.orfIdx < w.tr.orfs.length) { newProtein(); w.phase = 'scan'; w.t = 0; }
        else startPhase('riboOut');
      }
    }
  }

  // アミノ酸をのせたtRNAを1つ呼び寄せる（次のコドンを読むため）
  function callTRNA(a){
    const w = W, t = w.trnas[w.trnaNext++ % w.trnas.length];
    t.visible = true; t.scale = 1; t.state = 'in'; t.t = 0;
    t.quat.identity();
    t.pos.copy(w.ssu.pos).add(tmpA.set((Math.random() - 0.5) * 120, 40 + Math.random() * 40, 150));
    setCol(t.mesh, TRNA_AA, a.color); flush(t.mesh);
    w.trnaIn = t;
  }
  const TR_A = L(8, 14, 10), TR_OUT = L(-70, 50, 110);
  function stepTRNAs(dt, h){
    const w = W;
    w.trnas.forEach(t => {
      if (!t.visible) return;
      t.t += dt;
      t.update(h);
      if (t.state === 'in') {
        w.ssu.toWorld(TR_A, tmpA);
        t.pos.lerp(tmpA, 1 - Math.exp(-7 * dt));
        jiggle(t.pos, 0.4 * Math.sqrt(h));
      } else if (t.state === 'out') {
        w.ssu.toWorld(TR_OUT, tmpA);
        t.pos.lerp(tmpA, 1 - Math.exp(-1.6 * dt));
        jiggle(t.pos, 0.5 * Math.sqrt(h));
        t.scale = Math.max(0, 1 - t.t / 1.4);
        if (t.scale <= 0) { t.visible = false; t.state = 'idle'; }
      }
    });
  }

  function stepNTPs(dt, h){
    const w = W, sh = Math.sqrt(h);
    for (const nt of w.ntps) {
      if (nt.state !== 'free') { tumble(nt.q, 0.05 * sh); continue; }
      jiggle(nt.c, 1.5 * sh);
      tumble(nt.q, 0.18 * sh);
      for (const s of w.env.spheres) {
        tmpA.subVectors(nt.c, s.c);
        const d = tmpA.length();
        if (d < s.r + 4) nt.c.copy(s.c).addScaledVector(tmpA, (s.r + 4) / (d || 1e-6));
      }
      if (nt.c.x > w.env.x0 && nt.c.x < w.env.x1) {
        const ryz = Math.hypot(nt.c.y, nt.c.z) || 1e-6;
        if (ryz < 13) { nt.c.y *= 13 / ryz; nt.c.z *= 13 / ryz; }
      }
      if (nt.c.distanceTo(focus) > 125) respawnNTP(nt, focus, false);
    }
    for (let i = w.ppis.length - 1; i >= 0; i--) {
      const p = w.ppis[i];
      p.t += dt;
      p.c.addScaledVector(p.v, dt); p.v.multiplyScalar(Math.exp(-1.5 * dt));
      jiggle(p.c, 1.3 * sh); tumble(p.q, 0.2 * sh);
      if (p.t > 3) w.ppis.splice(i, 1);
    }
    for (const g of w.sigs) {
      jiggle(g.c, 1.2 * sh); tumble(g.q, 0.2 * sh);
      if (g.c.distanceTo(g.home) > 75) g.c.lerp(g.home, 0.05);
      const ryz = Math.hypot(g.c.y, g.c.z) || 1e-6;
      if (ryz < 12) { g.c.y *= 12 / ryz; g.c.z *= 12 / ryz; }
    }
  }

  // =====================================================================
  // 描画への反映
  // =====================================================================
  const NTP_OFF = { base: L(3.2, 0, 0), p1: L(-2.7, 1.0, 0), p2: L(-4.9, 2.0, 0.6), p3: L(-7.0, 2.9, 1.2) };
  function syncMeshes(){
    const w = W, h = w.helix, tr = w.tr;
    // DNA
    for (let i = 0; i < w.n; i++) {
      const o = i * 18, op = openAt(i), x = w.xOf(i);
      const ly = (k, ty, tz) => [h[o+k+1] + (ty - h[o+k+1]) * op, h[o+k+2] + (tz - h[o+k+2]) * op];
      let p;
      p = ly(0, 17.5, -3.5);  setInst(w.dP, 2*i,   x,     p[0], p[1], 3.0);
      p = ly(3, 14.5, -1.5);  setInst(w.dS, 2*i,   x+0.8, p[0], p[1], 2.4);
      p = ly(6, 11.5, 0.5);   setInst(w.dB, 2*i,   x+0.8, p[0], p[1], 2.45);
      p = ly(9, -13.5, -3.0); setInst(w.dP, 2*i+1, x,     p[0], p[1], 3.0);
      p = ly(12, -10.5, -0.5); setInst(w.dS, 2*i+1, x-0.8, p[0], p[1], 2.4);
      p = ly(15, -7.0, 2.0);  setInst(w.dB, 2*i+1, x-0.8, p[0], p[1], 2.45);
    }
    flush(w.dP); flush(w.dS); flush(w.dB);

    // RNA
    let nr = 0;
    if (w.rna) {
      const sc = w.rna.scale;
      w.rna.b.forEach((x, j) => {
        setInst(w.rP, j, x.p.x, x.p.y, x.p.z, 2.7 * sc);
        setInst(w.rB, j, x.q.x, x.q.y, x.q.z, 2.35 * sc);
        nr = j + 1;
      });
    }
    hideFrom(w.rP, nr); hideFrom(w.rB, nr);
    flush(w.rP); flush(w.rB);

    // NTP
    w.ntps.forEach((nt, i) => {
      setCol(w.nB, i, BASE_COL[nt.type]);
      tmpA.copy(NTP_OFF.base).applyQuaternion(nt.q).add(nt.c); setInst(w.nB, i, tmpA.x, tmpA.y, tmpA.z, 2.3);
      setInst(w.nS, i, nt.c.x, nt.c.y, nt.c.z, 2.0);
      ['p1','p2','p3'].forEach((k, m) => {
        tmpA.copy(NTP_OFF[k]).applyQuaternion(nt.q).add(nt.c); setInst(w.nP, i*3+m, tmpA.x, tmpA.y, tmpA.z, 1.55);
      });
    });
    let np = N_NTP * 3;
    w.ppis.forEach(p => {
      const s = 1.55 * Math.min(1, (3 - p.t) / 0.8);
      tmpA.set(1.25, 0, 0).applyQuaternion(p.q);
      setInst(w.nP, np++, p.c.x + tmpA.x, p.c.y + tmpA.y, p.c.z + tmpA.z, s);
      setInst(w.nP, np++, p.c.x - tmpA.x, p.c.y - tmpA.y, p.c.z - tmpA.z, s);
    });
    hideFrom(w.nP, np);
    flush(w.nB); flush(w.nS); flush(w.nP);

    let ns = 0;
    w.sigs.forEach(g => {
      tmpA.set(1.6, 0, 0).applyQuaternion(g.q);
      setInst(w.sigMesh, ns++, g.c.x + tmpA.x, g.c.y + tmpA.y, g.c.z + tmpA.z, 2.0);
      setInst(w.sigMesh, ns++, g.c.x - tmpA.x, g.c.y - tmpA.y, g.c.z - tmpA.z, 1.6);
    });
    hideFrom(w.sigMesh, ns); flush(w.sigMesh);

    // タンパク質
    let na = 0, ng = 0;
    const drawChain = c => {
      c.b.forEach(x => {
        if (x.glow) { if (ng < MAX_AA) { setInst(w.aaG, ng, x.p.x, x.p.y, x.p.z, 4.0 * c.scale); setCol(w.aaG, ng, x.color); ng++; } }
        else if (na < MAX_AA) { setInst(w.aa, na, x.p.x, x.p.y, x.p.z, 4.0 * c.scale); setCol(w.aa, na, x.color); na++; }
      });
    };
    if (w.prot) drawChain(w.prot);
    w.products.forEach(drawChain);
    hideFrom(w.aa, na); hideFrom(w.aaG, ng);
    flush(w.aa); flush(w.aaG);

    ['rnap','sigma','ssu','lsu','factor'].forEach(k => { if (w[k]) w[k].render(); });
    w.trnas.forEach(t => t.render());
  }

  // =====================================================================
  // カメラ・字幕・ミニマップ
  // =====================================================================
  function focusTarget(out){
    const w = W, ph = w.phase;
    if (ph === 'riboIn') return out.copy(w.ssu.pos).lerp(w.rna.b[Math.min(w.tr.orfs[0].start - w.tr.tStart, w.rna.b.length - 1)].p, 0.5);
    if (ph === 'lsuIn' || ph === 'translate' || ph === 'scan' || ph === 'riboOut') return out.copy(w.ssu.pos).add(tmpA.set(0, 18, 0));
    if (ph === 'fade') return w.products.length ? out.copy(w.products[w.products.length - 1].b[0].p) : out.copy(focus);
    if (ph === 'search' || ph === 'bounce') return out.set(0, 0, 0).lerp(w.rnap.pos, 0.35);
    // 転写中はポリメラーゼの少し後ろ(伸びていくmRNAの側)を中心に
    return out.copy(w.rnap.pos).add(tmpA.set(-14, 10, 0));
  }
  function phaseDist(){
    const ph = W.phase;
    if (ph === 'search' || ph === 'bounce') return 300;
    if (ph === 'elong' || ph === 'bind') return 215;
    if (ph === 'riboIn' || ph === 'fade') return 280;
    if (ph === 'translate' || ph === 'scan' || ph === 'lsuIn' || ph === 'riboOut') return 300;
    return 230;
  }

  function updateCamera(dt){
    focusTarget(camGoal);
    focus.lerp(camGoal, 1 - Math.exp(-1.6 * dt));
    let th, ph, dist;
    if (autoCam) {
      camT += dt;
      th = Math.PI / 2 - 0.45 + 0.42 * Math.sin(camT * 0.13);
      ph = 1.18 + 0.12 * Math.sin(camT * 0.09);
      userOrbit.dist += (phaseDist() - userOrbit.dist) * (1 - Math.exp(-0.8 * dt));
      dist = userOrbit.dist;
      userOrbit.theta = th; userOrbit.phi = ph;
    } else { th = userOrbit.theta; ph = userOrbit.phi; dist = userOrbit.dist; }
    camera.position.set(
      focus.x + dist * Math.sin(ph) * Math.cos(th),
      focus.y + dist * Math.cos(ph),
      focus.z + dist * Math.sin(ph) * Math.sin(th));
    camera.lookAt(focus);
  }

  function regionAt(bp){ return W.tr.regions.find(r => bp >= r.start && bp < r.end); }

  function blockedReason(){
    const k = W.circuit.promKind;
    if (k === 'prom-repressor') return 'リプレッサー(オレンジ)がプロモーターをふさいでいるため、RNAポリメラーゼが結合できません。';
    if (k === 'prom-activator') return 'アクチベーターが足りないため、このプロモーターにはRNAポリメラーゼがほとんど結合しません。';
    if (k === 'prom-signal-off') return 'シグナル分子(黄色)があるので、このプロモーターはOFFになっています。';
    if (k === 'prom-signal-on') return 'シグナル分子が足りないので、このプロモーターはOFFのままです。';
    return 'プロモーターの活性が低く、なかなか結合できません。';
  }

  function caption(){
    const w = W, tr = w.tr, ph = w.phase;
    const prom = '「' + tr.prom.label + '」';
    if (ph === 'search') return ['①', 'RNAポリメラーゼ(青い大きな分子)が、DNA上のプロモーター' + prom + 'を探しています。' +
      (w.fails && w.pBind < 0.5 ? '<br><small>' + blockedReason() + '</small>' : '')];
    if (ph === 'bounce') return ['✕', '結合できずに離れていきました。' + blockedReason()];
    if (ph === 'bind') return ['②', 'プロモーターに結合！ σ因子(黄土色)の助けでDNAの二重らせんがほどけ、<b>転写バブル</b>ができます。'];
    if (ph === 'elong') {
      const r = regionAt(w.s) || tr.prom;
      return ['③', '<b>転写</b>: 飛び回るヌクレオチド(NTP)のうち、鋳型DNAと対になるものだけが取り込まれ、mRNAが伸びていきます。' +
        '<br><small>いま読んでいる配列: 「' + (r.label || 'DNA') + '」 ／ DNAのT(黄)はRNAではU(ピンク)になります</small>'];
    }
    if (ph === 'term') return ['④', '<b>ターミネーター</b>: 転写されたRNAが<b>ヘアピン構造</b>をつくり、RNAポリメラーゼがDNAから外れます。'];
    if (ph === 'release') return tr.hasTerm
      ? ['④', 'mRNAが完成し、RNAポリメラーゼが離れていきます。DNAは元の二重らせんに戻ります。']
      : ['⚠️', 'ターミネーターがないため、DNAの端まで読み進めてしまいました（読み過ごし）。転写がきちんと終わらず、効率が下がります。'];
    if (ph === 'riboIn') return ['⑤', '<b>翻訳</b>の準備: リボソームの小サブユニット(うす紫)が、mRNAの開始コドン AUG を見つけて結合します。'];
    if (ph === 'lsuIn') return ['⑤', '大サブユニット(紫)が合体して、リボソームが完成します。'];
    if (ph === 'translate' || ph === 'scan') {
      const orf = w.prot && w.prot.orf;
      return ['⑥', 'リボソームがmRNAを<b>3塩基(コドン)ずつ</b>読みます。コドンに合うアミノ酸を<b>tRNA</b>(緑のL字)が運んできて、つながって「' + (orf ? orf.label : 'タンパク質') + '」になります。' +
        '<br><small>できた鎖はゆらぎながら折りたたまれ、立体構造になります' + (orf && orf.fusion ? '（2つのドメインがリンカーでつながった融合タンパク質）' : '') + '</small>'];
    }
    if (ph === 'riboOut') return ['⑦', '終止コドンに着くと、リボソームはばらばらになってmRNAから離れます。'];
    if (ph === 'fade') return ['⑦', tr.orfs.length
      ? '完成したタンパク質が細胞の中へ放出されました。mRNAはやがて分解され、この流れが何度もくり返されます。'
      : 'この回路には遺伝子(タンパク質コーディング領域)がないので、タンパク質は作られません。'];
    return ['', ''];
  }

  function buildMap(){
    if (!ui.map) return;
    const tr = W.tr;
    ui.map.innerHTML = tr.regions.filter(r => r.type !== 'pad').map(r =>
      '<span class="tx-seg" style="flex:' + (r.end - r.start) + ';background:#' + new THREE.Color(r.color).getHexString() + '" title="' + r.label + '">' +
        (r.short || '') + '</span>').join('') + '<i class="tx-marker" id="txMarker"></i>';
    ui.marker = ui.map.querySelector('#txMarker');
    const first = tr.regions[1].start, last = tr.regions[tr.regions.length - 2].end;
    ui.mapRange = [first, last];
  }

  let lastCap = '';
  function updateUI(){
    const w = W, c = caption(), key = c[0] + c[1];
    if (key !== lastCap && ui.cap) {
      lastCap = key;
      ui.cap.innerHTML = (c[0] ? '<span class="tx-step">' + c[0] + '</span>' : '') + '<span>' + c[1] + '</span>';
    }
    if (ui.marker) {
      const ph = w.phase;
      let bp = null, cls = 'rnap';
      if (ph === 'elong' || ph === 'bind' || ph === 'term') bp = w.sF;
      else if (ph === 'translate' || ph === 'scan' || ph === 'lsuIn') { bp = w.tr.tStart + w.kF; cls = 'ribo'; }
      if (bp == null) ui.marker.style.display = 'none';
      else {
        const [a, b] = ui.mapRange;
        ui.marker.style.display = 'block';
        ui.marker.className = 'tx-marker ' + cls;
        ui.marker.style.left = (Math.max(0, Math.min(1, (bp - a) / (b - a))) * 100).toFixed(2) + '%';
      }
    }
    if (ui.stats) {
      ui.stats.textContent = 'できたmRNA ' + w.mRNAs + ' 本 ／ タンパク質 ' + w.proteins + ' 個' +
        (w.fails ? ' ／ プロモーターへの結合: 成功 ' + w.binds + ' 回・失敗 ' + w.fails + ' 回' : '');
    }
  }

  // =====================================================================
  // ループ
  // =====================================================================
  function loop(ts){
    if (!running) return;
    // Labパネルを閉じていて見えないときは計算も描画も休む
    if (mountEl && mountEl.clientWidth === 0) { lastT = 0; raf = requestAnimationFrame(loop); return; }
    const realDt = lastT ? Math.min(0.05, (ts - lastT) / 1000) : 1 / 60;
    lastT = ts;
    const dt = realDt * speed;
    if (W && dt > 0) {
      const ns = Math.max(1, Math.ceil(dt / (1 / 45)));
      for (let k = 0; k < ns; k++) step(dt / ns);
    }
    if (W) {
      syncMeshes();
      updateCamera(realDt);
      backdrop.update(focus, realDt, camera);
      updateUI();
    }
    renderer.render(scene, camera);
    raf = requestAnimationFrame(loop);
  }

  function resize(){
    if (!renderer || !mountEl) return;
    const w = mountEl.clientWidth || 600, h = mountEl.clientHeight || 400;
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 1.6));
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  }

  function syncButtons(){
    if (!ui.root) return;
    ui.root.querySelectorAll('[data-speed]').forEach(b => b.classList.toggle('active', Number(b.dataset.speed) === speed));
    const cam = ui.root.querySelector('#txCam');
    if (cam) { cam.classList.toggle('active', autoCam); cam.textContent = autoCam ? '🎥 自動カメラ ON' : '🎥 自動カメラ OFF'; }
  }

  function buildDOM(container){
    container.innerHTML =
      '<div class="tx-root">' +
        '<div class="tx-shell" id="txShell">' +
          '<div class="tx-mount" id="txMount"></div>' +
          '<div class="tx-vignette"></div>' +
          '<div class="tx-map" id="txMap"></div>' +
          '<div class="tx-caption" id="txCap"></div>' +
          '<button class="tx-fs" id="txFs" title="全画面">⛶</button>' +
        '</div>' +
        '<div class="tx-bar">' +
          '<span class="tx-speed">速さ' +
            [0, 0.5, 1, 2, 4].map(s => '<button data-speed="' + s + '">' + (s === 0 ? '⏸' : '×' + s) + '</button>').join('') +
          '</span>' +
          '<button id="txCam">🎥 自動カメラ ON</button>' +
          '<button id="txRestart">⏮ 最初から</button>' +
        '</div>' +
        '<div class="tx-stats" id="txStats"></div>' +
      '</div>';
    ui = {
      root: container.querySelector('.tx-root'),
      shell: container.querySelector('#txShell'),
      map: container.querySelector('#txMap'),
      cap: container.querySelector('#txCap'),
      stats: container.querySelector('#txStats')
    };
    lastCap = '';
    ui.root.querySelectorAll('[data-speed]').forEach(b => b.addEventListener('click', () => { speed = Number(b.dataset.speed); syncButtons(); }));
    ui.root.querySelector('#txCam').addEventListener('click', () => { autoCam = !autoCam; syncButtons(); });
    ui.root.querySelector('#txRestart').addEventListener('click', () => { if (W) buildWorld(W.circuit); });
    ui.root.querySelector('#txFs').addEventListener('click', () => {
      const el = ui.shell;
      if (document.fullscreenElement) document.exitFullscreen();
      else if (el.requestFullscreen) el.requestFullscreen();
    });
    syncButtons();
    return container.querySelector('#txMount');
  }

  // circuit = { blocks:[part...], promKind, k_tx, signal, label }
  function mount(container, circuit){
    if (typeof THREE === 'undefined') return false;
    if (!renderer) initRenderer();
    host = container;
    mountEl = buildDOM(container);
    mountEl.appendChild(renderer.domElement);
    if (ro) ro.disconnect();
    if (window.ResizeObserver) { ro = new ResizeObserver(resize); ro.observe(mountEl); }
    buildWorld(circuit);
    resize();
    if (!running) { running = true; lastT = 0; raf = requestAnimationFrame(loop); }
    return true;
  }

  function unmount(){
    running = false;
    if (raf) cancelAnimationFrame(raf);
    raf = null;
    if (ro) { ro.disconnect(); ro = null; }
  }

  // _advance はテスト用: 描画を待たずにシミュレーションを sec 秒進める
  window.TxAnim = { mount, unmount, _debug: () => W,
    _advance: sec => {
      for (let t = 0; t < sec && W; t += 1 / 45) step(1 / 45);
      if (W) { focusTarget(focus); userOrbit.dist = phaseDist(); }
    } };
})();
