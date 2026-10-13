/*
 * 酵素のドメイン運動と触媒サイクル（粗視化分子動力学）
 *
 * アミノ酸1残基 = 1粒子(Cα)に簡略化した粗視化モデルで、酵素分子の
 * 「可動部分」を計算します。
 *
 *   骨格      : 結合長・1-3・1-4 距離の調和バネ（ドメイン内のみ）
 *   立体構造  : 天然コンタクト(Go模型・モース型ポテンシャル)
 *   排除体積  : 全粒子間の斥力
 *   熱運動    : ランジュバン熱浴（摩擦＋熱ゆらぎ）
 *
 * 可動部分の作り方:
 *   - リンカー部には 1-3/1-4 拘束も天然コンタクトも入れません。
 *     → 2つのドメインがリンカー上で自由に振れる（ヒンジ運動）
 *   - 触媒ドメインは2ローブ構成にし、裂け目の「底」付近だけを
 *     天然コンタクトで留めます。
 *     → 活性部位が熱ゆらぎで開閉する（呼吸運動）
 *
 * 触媒サイクル:
 *   基質(アルカン鎖)が拡散 → 活性部位に捕捉 → 一定時間保持 → 切断 →
 *   生成物が離れる → 新しい基質が供給される
 *
 * 実在のアミノ酸配列ではなく、簡略化した仮想の構造です。
 */
(function(){
  'use strict';

  // --- 力場パラメータ（既存の3Dタンパク質MDと同じスケール）--------------
  const K_BOND = 80, K13 = 25, K14 = 8;
  const SIGMA_REP = 2.8, EPS_REP = 0.8;
  const EPS_NATIVE = 4.5, MORSE_A = 1.3;
  const F_CAP = 60, MAX_V = 25;
  const GAMMA = 3.0, DT = 0.006;
  // 基質は低分子なので摩擦が小さく、拡散が速い（D = kT/γ）
  const GAMMA_SUB = 0.5;
  const CUTOFF = 7.2;
  const BOND = 3.8;

  // --- 触媒サイクル ------------------------------------------------------
  const SUB_LEN = 8;              // 基質アルカン鎖の粒子数
  const R_CAPTURE = 18;           // 活性部位が基質を引き寄せ始める距離
  const R_BOUND = 6.5;            // 結合とみなす距離
  const EPS_BIND = 6.0;           // 捕捉ポテンシャルの深さ
  const DWELL = 60;               // 切断までの保持ステップ数
  const BOX = 34;                 // 基質が拡散する領域の半径
  const SPAWN_R = 0.55;           // 供給位置（BOX に対する比）

  const clamp = (v, a, b) => Math.max(a, Math.min(b, v));

  // 構造生成用のシード付き乱数。同じシードなら必ず同じ立体構造になるので、
  // 「安定性だけを変えて比べる」ことができます（熱ゆらぎ側は Math.random）。
  let _seed = 1;
  function srand(){ _seed = (_seed * 1664525 + 1013904223) >>> 0; return _seed / 4294967296; }
  function setSeed(s){ _seed = (s >>> 0) || 1; }
  function gauss(){
    let u = 0, v = 0;
    while (u === 0) u = Math.random();
    while (v === 0) v = Math.random();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }
  const dist = (a, b) => Math.hypot(b.x - a.x, b.y - a.y, b.z - a.z);

  // =====================================================================
  // 構造の生成
  // =====================================================================
  // 球内に収まるコンパクトな自己回避鎖（毛糸玉状のドメイン）。
  // 半径は粒子数から決めます。最小間隔 0.85·BOND の粒子 n 個を自己回避ウォークで
  // 詰めるには充填率 2 割程度が上限なので、それを満たす半径を使います。
  const MIN_SEP = 0.85;
  function globuleRadius(n){ return 0.85 * BOND * Math.cbrt(n); }

  function globule(n, center, seedDir){
    const R = globuleRadius(n);
    const pts = [];
    let cur = { x: center.x, y: center.y, z: center.z };
    let dir = seedDir || { x: 1, y: 0, z: 0 };
    const clampToSphere = p => {
      const dx=p.x-center.x, dy=p.y-center.y, dz=p.z-center.z;
      const r=Math.hypot(dx,dy,dz);
      if (r <= R) return p;
      const k=R/r;
      return { x:center.x+dx*k, y:center.y+dy*k, z:center.z+dz*k };
    };
    for (let i = 0; i < n; i++) {
      if (i === 0) { pts.push({ x: cur.x, y: cur.y, z: cur.z }); continue; }
      let placed = null;
      for (let attempt = 0; attempt < 220; attempt++) {
        const nd = {
          x: dir.x + (srand() - 0.5) * 1.8,
          y: dir.y + (srand() - 0.5) * 1.8,
          z: dir.z + (srand() - 0.5) * 1.8
        };
        const toC = { x: center.x - cur.x, y: center.y - cur.y, z: center.z - cur.z };
        const rad = Math.hypot(toC.x, toC.y, toC.z);
        if (rad > R * 0.6) {
          const w = (rad - R * 0.6) / R * 2.5;
          nd.x += toC.x / rad * w; nd.y += toC.y / rad * w; nd.z += toC.z / rad * w;
        }
        const nn = Math.hypot(nd.x, nd.y, nd.z) || 1;
        nd.x /= nn; nd.y /= nn; nd.z /= nn;
        const cand = { x: cur.x + nd.x * BOND, y: cur.y + nd.y * BOND, z: cur.z + nd.z * BOND };
        if (dist(cand, center) > R) continue;
        let ok = true;
        for (let k = 0; k < pts.length - 1; k++) {
          if (dist(pts[k], cand) < BOND * MIN_SEP) { ok = false; break; }
        }
        if (ok) { placed = cand; dir = nd; break; }
      }
      if (!placed) {
        // 行き詰まったら球内に引き戻す（球外へ逃げないことを保証する）
        placed = clampToSphere({ x: cur.x + dir.x * BOND, y: cur.y + dir.y * BOND, z: cur.z + dir.z * BOND });
      }
      pts.push(placed);
      cur = placed;
    }
    return pts;
  }

  // リンカー: 2点を結ぶ緩い鎖（拘束を入れないので自由に曲がる）
  function linkerChain(from, to, n){
    const pts = [];
    for (let i = 1; i <= n; i++) {
      const t = i / (n + 1);
      pts.push({
        x: from.x + (to.x - from.x) * t + (srand() - 0.5) * 2,
        y: from.y + (to.y - from.y) * t + (srand() - 0.5) * 2,
        z: from.z + (to.z - from.z) * t + (srand() - 0.5) * 2
      });
    }
    return pts;
  }

  // 酵素の設計図を組み立てる。
  // twoDomain=true なら「触媒ドメイン＋リンカー＋安定化ドメイン」（融合酵素）
  function buildEnzyme(opts){
    setSeed(opts.seed != null ? opts.seed : 12345);
    const twoDomain = !!opts.twoDomain;
    const hasCleft = opts.hasCleft !== false;
    const pos = [], domain = [], name = [];
    const push = (p, d) => { pos.push(p); domain.push(d); };

    // --- 触媒ドメイン: 2ローブ + ヒンジ。ローブの間が活性部位の裂け目 ---
    const N_LOBE = 15, N_STAB = 24;
    const lobeR = globuleRadius(N_LOBE), sep = lobeR * 0.75;
    const lobeA = globule(N_LOBE, { x:-sep, y:0, z:0 }, { x:0, y:-1, z:0 });
    lobeA.forEach(p => push(p, 'lobeA'));
    const idxHingeStart = pos.length;
    // ヒンジは裂け目の「底」(y<0)を通す
    const hinge = linkerChain(lobeA[lobeA.length-1], { x: sep, y: -lobeR * 0.75, z: 0 }, 4)
      .map(p => ({ x: p.x, y: Math.min(p.y, -lobeR * 0.55), z: p.z }));
    hinge.forEach(p => push(p, 'hinge'));
    const lobeB = globule(N_LOBE, { x: sep, y:0, z:0 }, { x:0, y:1, z:0 });
    lobeB.forEach(p => push(p, 'lobeB'));
    const catEnd = pos.length - 1;

    let linkStart = -1, linkEnd = -1, stabStart = -1;
    if (twoDomain) {
      const stabCenter = { x: 1, y: lobeR * 2.5, z: 4 };
      linkStart = pos.length;
      linkerChain(pos[catEnd], stabCenter, 6).forEach(p => push(p, 'linker'));
      linkEnd = pos.length - 1;
      stabStart = pos.length;
      globule(N_STAB, stabCenter, { x:0, y:1, z:0 }).forEach(p => push(p, 'stab'));
    }

    const N = pos.length;

    // --- 活性部位: 裂け目の上部(y>0)に面したローブの粒子 ---
    const active = [];
    for (let i = 0; i < N; i++) {
      if (domain[i] !== 'lobeA' && domain[i] !== 'lobeB') continue;
      if (pos[i].y > lobeR * 0.15 && Math.abs(pos[i].x) < sep + 2.5) active.push(i);
    }

    // --- トポロジー ---
    const sameDomainGroup = i =>
      (domain[i] === 'lobeA' || domain[i] === 'lobeB' || domain[i] === 'hinge') ? 'cat' : domain[i];

    const bondLen = [], d13 = [], d14 = [], contacts = [];
    for (let i = 0; i < N - 1; i++) bondLen.push(dist(pos[i], pos[i+1]));
    for (let i = 0; i < N - 2; i++) {
      // リンカーを跨ぐ角度拘束は入れない → ヒンジが自由に曲がる
      const free = domain[i] === 'linker' || domain[i+1] === 'linker' || domain[i+2] === 'linker';
      d13.push(free ? -1 : dist(pos[i], pos[i+2]));
    }
    for (let i = 0; i < N - 3; i++) {
      const free = domain[i] === 'linker' || domain[i+1] === 'linker' ||
                   domain[i+2] === 'linker' || domain[i+3] === 'linker';
      d14.push(free ? -1 : dist(pos[i], pos[i+3]));
    }

    const stabScale = opts.stability != null ? opts.stability : 1;
    // ローブ粒子のy座標の中央値。これより下を「裂け目の底」とみなす。
    const lobeYs = [];
    for (let i = 0; i < N; i++) if (domain[i]==='lobeA'||domain[i]==='lobeB') lobeYs.push(pos[i].y);
    lobeYs.sort((a,b)=>a-b);
    const yMid = lobeYs.length ? lobeYs[Math.floor(lobeYs.length/2)] : 0;
    for (let i = 0; i < N; i++) {
      for (let j = i + 4; j < N; j++) {
        const d = dist(pos[i], pos[j]);
        if (d >= CUTOFF) continue;
        const gi = sameDomainGroup(i), gj = sameDomainGroup(j);
        if (gi === 'linker' || gj === 'linker') continue;   // リンカーは拘束しない
        if (gi !== gj) continue;                            // ドメイン間も拘束しない（＝ヒンジ）
        const cross = (domain[i] === 'lobeA' && domain[j] === 'lobeB') ||
                      (domain[i] === 'lobeB' && domain[j] === 'lobeA');
        if (cross) continue;                                // 下で明示的に張る
        const eps = (gi === 'stab') ? 1.25 : 1;             // 安定化ドメインは密に留める
        contacts.push({ i, j, d0: d, eps: eps * stabScale });
      }
    }


    // --- 裂け目のアンカー ---
    // ローブ間は距離カットオフでは安定して拾えない（生成のばらつきが大きい）。
    // 裂け目の「底」側にある最近接ペアだけを明示的に留め、上部は開放したままにする。
    // これで活性部位は形を保ちつつ熱ゆらぎで開閉できる。
    const crossPairs = [];
    for (let i = 0; i < N; i++) {
      if (domain[i] !== 'lobeA') continue;
      for (let j = 0; j < N; j++) {
        if (domain[j] !== 'lobeB') continue;
        if (pos[i].y >= yMid || pos[j].y >= yMid) continue;   // 上部は留めない
        crossPairs.push({ i, j, d: dist(pos[i], pos[j]) });
      }
    }
    crossPairs.sort((a, b) => a.d - b.d);
    const nAnchor = hasCleft ? 6 : 12;
    crossPairs.slice(0, nAnchor).forEach(p => {
      contacts.push({ i: p.i, j: p.j, d0: p.d, eps: 0.8 * stabScale });
    });

    return {
      N, pos, domain, active, bondLen, d13, d14, contacts,
      twoDomain, catEnd, linkStart, linkEnd, stabStart,
      nativeRg: rgOf(pos)
    };
  }

  function rgOf(p){
    const n = p.length;
    let cx=0, cy=0, cz=0;
    p.forEach(q => { cx+=q.x; cy+=q.y; cz+=q.z; });
    cx/=n; cy/=n; cz/=n;
    let s = 0;
    p.forEach(q => { s += (q.x-cx)**2 + (q.y-cy)**2 + (q.z-cz)**2; });
    return Math.sqrt(s/n);
  }
  function centroid(pos, pick){
    let cx=0, cy=0, cz=0, n=0;
    for (let i = 0; i < pos.length; i++) {
      if (pick && !pick(i)) continue;
      cx+=pos[i].x; cy+=pos[i].y; cz+=pos[i].z; n++;
    }
    return n ? { x:cx/n, y:cy/n, z:cz/n } : { x:0, y:0, z:0 };
  }

  // =====================================================================
  // シミュレーション
  // =====================================================================
  function Engine(spec){
    this.spec = spec;
    this.pos = spec.pos.map(p => ({ x:p.x, y:p.y, z:p.z }));
    this.vel = this.pos.map(() => ({ x:0, y:0, z:0 }));
    this.kT = 1.2;
    this.substrates = [];
    this.turnovers = 0;
    this.steps = 0;
    this.lastCleave = -1e9;
    for (let i = 0; i < 5; i++) this.spawnSubstrate();
  }

  Engine.prototype.spawnSubstrate = function(){
    // 箱の縁からランダムな向きで供給する
    const th = Math.random()*Math.PI*2, ph = Math.acos(Math.random()*2-1);
    const R = BOX*SPAWN_R;
    const c = {
      x: R*Math.sin(ph)*Math.cos(th),
      y: R*Math.cos(ph) + 8,
      z: R*Math.sin(ph)*Math.sin(th)
    };
    const d = { x:Math.random()-0.5, y:Math.random()-0.5, z:Math.random()-0.5 };
    const nn = Math.hypot(d.x,d.y,d.z)||1; d.x/=nn; d.y/=nn; d.z/=nn;
    const beads = [];
    for (let i = 0; i < SUB_LEN; i++) {
      beads.push({ x:c.x+d.x*BOND*i, y:c.y+d.y*BOND*i, z:c.z+d.z*BOND*i,
                   vx:0, vy:0, vz:0 });
    }
    this.substrates.push({ beads, state:'free', dwell:0, cut:-1, age:0 });
  };

  // 活性部位の位置は、2ローブの中点からヒンジの反対方向へずらした点とします。
  // 粒子の重心そのものではなく方向ベクトルで定義するので、分子が回転しても
  // 裂け目の開口部を正しく指し続けます。
  Engine.prototype.activeSite = function(){
    const sp = this.spec, pos = this.pos;
    const a = centroid(pos, i => sp.domain[i] === 'lobeA');
    const b = centroid(pos, i => sp.domain[i] === 'lobeB');
    const h = centroid(pos, i => sp.domain[i] === 'hinge');
    const mid = { x:(a.x+b.x)/2, y:(a.y+b.y)/2, z:(a.z+b.z)/2 };
    const up = { x:mid.x-h.x, y:mid.y-h.y, z:mid.z-h.z };
    const n = Math.hypot(up.x, up.y, up.z) || 1;
    const off = 2.5;
    return { x: mid.x + up.x/n*off, y: mid.y + up.y/n*off, z: mid.z + up.z/n*off };
  };

  // 活性部位の開き（2ローブの重心間距離）
  Engine.prototype.cleftWidth = function(){
    const a = centroid(this.pos, i => this.spec.domain[i] === 'lobeA');
    const b = centroid(this.pos, i => this.spec.domain[i] === 'lobeB');
    return dist(a, b);
  };

  // ドメイン間の角度（ヒンジの振れ）
  Engine.prototype.hingeAngle = function(){
    if (!this.spec.twoDomain) return null;
    const cat = centroid(this.pos, i => this.spec.domain[i] === 'lobeA' || this.spec.domain[i] === 'lobeB');
    const stab = centroid(this.pos, i => this.spec.domain[i] === 'stab');
    const link = centroid(this.pos, i => this.spec.domain[i] === 'linker');
    const v1 = { x:cat.x-link.x, y:cat.y-link.y, z:cat.z-link.z };
    const v2 = { x:stab.x-link.x, y:stab.y-link.y, z:stab.z-link.z };
    const n1 = Math.hypot(v1.x,v1.y,v1.z)||1, n2 = Math.hypot(v2.x,v2.y,v2.z)||1;
    const c = clamp((v1.x*v2.x + v1.y*v2.y + v1.z*v2.z)/(n1*n2), -1, 1);
    return Math.acos(c) * 180 / Math.PI;
  };

  // 天然コンタクト形成率（構造がどれだけ保たれているか）
  Engine.prototype.Q = function(){
    const c = this.spec.contacts;
    if (!c.length) return 1;
    let f = 0;
    c.forEach(k => { if (dist(this.pos[k.i], this.pos[k.j]) < k.d0 * 1.2) f++; });
    return f / c.length;
  };

  const cap = f => clamp(f, -F_CAP, F_CAP);

  Engine.prototype.step = function(){
    const sp = this.spec, N = sp.N, pos = this.pos, vel = this.vel;
    const F = pos.map(() => ({ x:0, y:0, z:0 }));

    const spring = (i, j, r0, k) => {
      if (r0 < 0) return;
      const dx = pos[j].x-pos[i].x, dy = pos[j].y-pos[i].y, dz = pos[j].z-pos[i].z;
      const r = Math.hypot(dx,dy,dz) || 1e-6;
      const m = cap(k*(r-r0)/r);
      F[i].x += m*dx; F[i].y += m*dy; F[i].z += m*dz;
      F[j].x -= m*dx; F[j].y -= m*dy; F[j].z -= m*dz;
    };
    for (let i = 0; i < N-1; i++) spring(i, i+1, sp.bondLen[i], K_BOND);
    for (let i = 0; i < N-2; i++) spring(i, i+2, sp.d13[i], K13);
    for (let i = 0; i < N-3; i++) spring(i, i+3, sp.d14[i], K14);

    // 天然コンタクト（モース型: 熱で切れる）
    sp.contacts.forEach(c => {
      const a = pos[c.i], b = pos[c.j];
      const dx = b.x-a.x, dy = b.y-a.y, dz = b.z-a.z;
      const r = Math.hypot(dx,dy,dz) || 1e-6;
      const ex = Math.exp(-MORSE_A*(r-c.d0));
      // U = eps(1-ex)^2,  dU/dr = 2·eps·a·ex·(1-ex),  F_i = (dU/dr)·d̂
      const m = cap((2*EPS_NATIVE*c.eps*MORSE_A*ex*(1-ex))/r);
      F[c.i].x += m*dx; F[c.i].y += m*dy; F[c.i].z += m*dz;
      F[c.j].x -= m*dx; F[c.j].y -= m*dy; F[c.j].z -= m*dz;
    });

    // 排除体積
    for (let i = 0; i < N; i++) {
      for (let j = i+3; j < N; j++) {
        const dx = pos[j].x-pos[i].x, dy = pos[j].y-pos[i].y, dz = pos[j].z-pos[i].z;
        const r2 = dx*dx+dy*dy+dz*dz;
        if (r2 >= SIGMA_REP*SIGMA_REP) continue;
        const r = Math.sqrt(r2) || 1e-6;
        const sr6 = Math.pow(SIGMA_REP/r, 6);
        const m = cap(12*EPS_REP*sr6*sr6/r/r);
        F[i].x -= m*dx; F[i].y -= m*dy; F[i].z -= m*dz;
        F[j].x += m*dx; F[j].y += m*dy; F[j].z += m*dz;
      }
    }

    const noise = Math.sqrt(2*GAMMA*this.kT*DT);
    for (let i = 0; i < N; i++) {
      ['x','y','z'].forEach(ax => {
        vel[i][ax] += (F[i][ax] - GAMMA*vel[i][ax])*DT + noise*gauss();
        vel[i][ax] = clamp(vel[i][ax], -MAX_V, MAX_V);
        pos[i][ax] += vel[i][ax]*DT;
      });
    }

    this.stepSubstrates();
    this.steps++;
  };

  Engine.prototype.stepSubstrates = function(){
    const site = this.activeSite();
    const noise = Math.sqrt(2*GAMMA_SUB*this.kT*DT);
    const prot = this.pos, sp = this.spec;

    this.substrates.forEach(sub => {
      const b = sub.beads;
      sub.age++;
      const mid = b[Math.floor(b.length/2)];
      const dSite = Math.hypot(mid.x-site.x, mid.y-site.y, mid.z-site.z);

      // 結合ポテンシャル（活性部位に向かう引力）
      const attract = sub.state === 'free' && dSite < R_CAPTURE;
      for (let i = 0; i < b.length; i++) {
        let fx = 0, fy = 0, fz = 0;
        // 鎖の結合
        if (i > 0) {
          const dx=b[i-1].x-b[i].x, dy=b[i-1].y-b[i].y, dz=b[i-1].z-b[i].z;
          const r=Math.hypot(dx,dy,dz)||1e-6;
          if (!(sub.cut >= 0 && i === sub.cut)) {
            const m=cap(K_BOND*0.5*(r-BOND)/r); fx+=m*dx; fy+=m*dy; fz+=m*dz;
          }
        }
        if (i < b.length-1) {
          const dx=b[i+1].x-b[i].x, dy=b[i+1].y-b[i].y, dz=b[i+1].z-b[i].z;
          const r=Math.hypot(dx,dy,dz)||1e-6;
          if (!(sub.cut >= 0 && i+1 === sub.cut)) {
            const m=cap(K_BOND*0.5*(r-BOND)/r); fx+=m*dx; fy+=m*dy; fz+=m*dz;
          }
        }
        // タンパク質との排除体積
        for (let k = 0; k < sp.N; k++) {
          const dx=prot[k].x-b[i].x, dy=prot[k].y-b[i].y, dz=prot[k].z-b[i].z;
          const r2=dx*dx+dy*dy+dz*dz;
          if (r2 >= SIGMA_REP*SIGMA_REP) continue;
          const r=Math.sqrt(r2)||1e-6;
          const sr6=Math.pow(SIGMA_REP/r,6);
          const m=cap(12*EPS_REP*sr6*sr6/r/r);
          fx-=m*dx; fy-=m*dy; fz-=m*dz;
        }
        if (attract) {
          const dx=site.x-b[i].x, dy=site.y-b[i].y, dz=site.z-b[i].z;
          const r=Math.hypot(dx,dy,dz)||1e-6;
          const m=EPS_BIND*Math.exp(-r/R_CAPTURE);
          fx+=m*dx/r*4; fy+=m*dy/r*4; fz+=m*dz/r*4;
        }
        // 領域外に出たら戻す
        const rr = Math.hypot(b[i].x, b[i].y-8, b[i].z);
        if (rr > BOX) {
          const m = -(rr-BOX)*0.6/rr;
          fx += m*b[i].x; fy += m*(b[i].y-8); fz += m*b[i].z;
        }
        b[i].vx += (fx - GAMMA_SUB*b[i].vx)*DT + noise*gauss();
        b[i].vy += (fy - GAMMA_SUB*b[i].vy)*DT + noise*gauss();
        b[i].vz += (fz - GAMMA_SUB*b[i].vz)*DT + noise*gauss();
        b[i].vx = clamp(b[i].vx,-MAX_V,MAX_V);
        b[i].vy = clamp(b[i].vy,-MAX_V,MAX_V);
        b[i].vz = clamp(b[i].vz,-MAX_V,MAX_V);
        b[i].x += b[i].vx*DT; b[i].y += b[i].vy*DT; b[i].z += b[i].vz*DT;
      }

      // 触媒サイクル
      if (sub.state === 'free') {
        if (dSite < R_BOUND) {
          sub.dwell++;
          if (sub.dwell > DWELL) {
            sub.state = 'product';
            sub.cut = Math.floor(b.length/2);   // 鎖の中央で切断
            this.turnovers++;
            this.lastCleave = this.steps;
          }
        } else {
          sub.dwell = Math.max(0, sub.dwell - 2);
        }
      }
    });

    // 十分離れた生成物は取り除き、新しい基質を供給する
    for (let i = this.substrates.length - 1; i >= 0; i--) {
      const s = this.substrates[i];
      if (s.state !== 'product') continue;
      const mid = s.beads[0];
      if (Math.hypot(mid.x, mid.y-8, mid.z) > BOX*0.75 && s.age > 150) {
        this.substrates.splice(i, 1);
        this.spawnSubstrate();
      }
    }
    while (this.substrates.length < 5) this.spawnSubstrate();
  };

  Engine.prototype.reset = function(){
    this.pos = this.spec.pos.map(p => ({
      x:p.x + (Math.random()-0.5)*1.2,
      y:p.y + (Math.random()-0.5)*1.2,
      z:p.z + (Math.random()-0.5)*1.2
    }));
    this.vel = this.pos.map(() => ({ x:0, y:0, z:0 }));
    this.substrates = [];
    this.turnovers = 0; this.steps = 0;
    for (let i = 0; i < 5; i++) this.spawnSubstrate();
  };

  window.EnzymeMD = { buildEnzyme: buildEnzyme, Engine: Engine, BOX: BOX, SUB_LEN: SUB_LEN, rgOf: rgOf };
})();

/* =====================================================================
 * 描画（Three.js）
 * ===================================================================== */
(function(){
  'use strict';
  if (!window.EnzymeMD) return;
  const { buildEnzyme, Engine, BOX } = window.EnzymeMD;

  const COLOR = {
    lobeA:  0xf97316, lobeB: 0xfb923c, hinge: 0xb45309,
    linker: 0x94a3b8, stab:  0x14b8a6
  };
  const C_ACTIVE = 0xfde047, C_SUB = 0x3f3f46, C_PROD = 0xa1a1aa;
  // 粒子を結合長(3.8)より大きめにして重ね、でこぼこした分子表面のように見せる
  const BEAD_R = 2.35, BOND_R = 0.9;
  const ML = window.MolLook;
  const mat = (hex, o) => ML ? ML.material(hex, o) :
    new THREE.MeshStandardMaterial(Object.assign({ color: hex, roughness: 0.5 }, o && o.emissive != null ? { emissive: o.emissive } : {}));

  let host=null, renderer=null, scene=null, camera=null, raf=null, running=false;
  let eng=null, spec=null;
  let beadMeshes={}, bondMesh=null, subMesh=null, siteMesh=null;
  let orbit={theta:0.9, phi:1.15, dist:95}, dragging=false, lastPtr=null;
  let viewTarget={x:0,y:0,z:0};
  let substeps=3, flash=0, tRate=[], onStats=null, backdrop=null, lastFrame=0;

  const M4 = () => new THREE.Matrix4();

  function makeScene(){
    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(45, 1, 0.5, 900);
    if (ML) {
      ML.setupScene(scene, 0.0075);
      backdrop = new ML.Backdrop(scene, 30, 150);
    } else {
      scene.background = new THREE.Color(0x0a1210);
      scene.add(new THREE.AmbientLight(0xffffff, 0.55));
      const l1=new THREE.DirectionalLight(0xffffff,0.95); l1.position.set(22,30,26); scene.add(l1);
    }
  }

  function disposeMesh(m){ if(!m) return; scene.remove(m); m.geometry.dispose(); m.material.dispose(); }

  function buildMeshes(){
    Object.values(beadMeshes).forEach(disposeMesh); beadMeshes={};
    disposeMesh(bondMesh); disposeMesh(subMesh); disposeMesh(siteMesh);

    const groups={};
    spec.domain.forEach((d,i)=>{ (groups[d]=groups[d]||[]).push(i); });
    const isActive = new Set(spec.active);
    Object.keys(groups).forEach(key=>{
      const idx = groups[key].filter(i=>!isActive.has(i));
      if(!idx.length) return;
      const geo=new THREE.IcosahedronGeometry(BEAD_R,2);
      const m=new THREE.InstancedMesh(geo,mat(COLOR[key]||0x999999,{roughness:0.6,rimStr:0.6}),idx.length);
      m.frustumCulled=false; m.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
      m.userData.idx=idx; scene.add(m); beadMeshes[key]=m;
    });
    // 活性部位は発光させて目立たせる
    if(spec.active.length){
      const geo=new THREE.IcosahedronGeometry(BEAD_R*1.08,2);
      siteMesh=new THREE.InstancedMesh(geo,mat(C_ACTIVE,{emissive:0x8a6d00,roughness:0.35,rim:0xfff3b0,rimStr:0.8}),spec.active.length);
      siteMesh.frustumCulled=false; siteMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
      siteMesh.userData.idx=spec.active.slice(); scene.add(siteMesh);
    }
    // 骨格
    const bgeo=new THREE.CylinderGeometry(BOND_R,BOND_R,1,7,1);
    bondMesh=new THREE.InstancedMesh(bgeo,mat(0xcbd5e1,{roughness:0.6,rimStr:0.3}),spec.N-1);
    bondMesh.frustumCulled=false; bondMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    scene.add(bondMesh);
    // 基質
    const sgeo=new THREE.IcosahedronGeometry(BEAD_R*0.78,1);
    subMesh=new THREE.InstancedMesh(sgeo,mat(0xffffff,{vertexColor:true,roughness:0.5,rim:0xffd9a0,rimStr:0.7}),12*window.EnzymeMD.SUB_LEN);
    subMesh.frustumCulled=false;
    subMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    subMesh.instanceColor=new THREE.InstancedBufferAttribute(new Float32Array(12*window.EnzymeMD.SUB_LEN*3),3);
    scene.add(subMesh);
  }

  const _m=new THREE.Matrix4(), _q=new THREE.Quaternion(), _v=new THREE.Vector3(),
        _up=new THREE.Vector3(0,1,0), _s=new THREE.Vector3(1,1,1), _c=new THREE.Color(), _vt=new THREE.Vector3();

  function sync(){
    const pos=eng.pos;
    Object.values(beadMeshes).forEach(mesh=>{
      mesh.userData.idx.forEach((pi,n)=>{
        _m.makeTranslation(pos[pi].x,pos[pi].y,pos[pi].z);
        mesh.setMatrixAt(n,_m);
      });
      mesh.instanceMatrix.needsUpdate=true;
    });
    if(siteMesh){
      siteMesh.userData.idx.forEach((pi,n)=>{
        _m.makeTranslation(pos[pi].x,pos[pi].y,pos[pi].z);
        siteMesh.setMatrixAt(n,_m);
      });
      siteMesh.instanceMatrix.needsUpdate=true;
      siteMesh.material.emissiveIntensity = 1 + flash*2.5;
    }
    for(let i=0;i<spec.N-1;i++){
      const a=pos[i], b=pos[i+1];
      _v.set(b.x-a.x,b.y-a.y,b.z-a.z);
      const len=_v.length()||1e-6;
      _q.setFromUnitVectors(_up,_v.clone().normalize());
      _s.set(1,len,1);
      _m.compose(new THREE.Vector3((a.x+b.x)/2,(a.y+b.y)/2,(a.z+b.z)/2),_q,_s);
      bondMesh.setMatrixAt(i,_m);
    }
    bondMesh.instanceMatrix.needsUpdate=true;

    let n=0;
    eng.substrates.forEach(sub=>{
      const col = sub.state==='product' ? C_PROD : C_SUB;
      sub.beads.forEach(b=>{
        if(n>=subMesh.count) return;
        _m.makeTranslation(b.x,b.y,b.z);
        subMesh.setMatrixAt(n,_m);
        _c.setHex(col); if (ML) _c.convertSRGBToLinear(); subMesh.setColorAt(n,_c);
        n++;
      });
    });
    for(let k=n;k<subMesh.count;k++){ _m.makeTranslation(0,-9999,0); subMesh.setMatrixAt(k,_m); }
    subMesh.instanceMatrix.needsUpdate=true;
    if(subMesh.instanceColor) subMesh.instanceColor.needsUpdate=true;
  }

  function loop(){
    if(!running) return;
    const before=eng.turnovers;
    for(let i=0;i<substeps;i++) eng.step();
    if(eng.turnovers>before){ flash=1; tRate.push(performance.now()); }
    flash*=0.90;
    while(tRate.length && performance.now()-tRate[0]>10000) tRate.shift();
    sync();
    const cx=orbit.dist*Math.sin(orbit.phi)*Math.cos(orbit.theta);
    const cy=orbit.dist*Math.cos(orbit.phi);
    const cz=orbit.dist*Math.sin(orbit.phi)*Math.sin(orbit.theta);
    camera.position.set(viewTarget.x+cx, viewTarget.y+cy, viewTarget.z+cz);
    camera.lookAt(viewTarget.x, viewTarget.y, viewTarget.z);
    const now=performance.now(), fdt=lastFrame?Math.min(0.05,(now-lastFrame)/1000):1/60; lastFrame=now;
    if(backdrop) backdrop.update(_vt.set(viewTarget.x,viewTarget.y,viewTarget.z), fdt, camera);
    renderer.render(scene,camera);
    if(onStats) onStats({
      cleft: eng.cleftWidth(), hinge: eng.hingeAngle(), Q: eng.Q(),
      turnovers: eng.turnovers, rate: tRate.length*6, kT: eng.kT
    });
    raf=requestAnimationFrame(loop);
  }

  function attach(dom){
    dom.addEventListener('pointerdown',e=>{dragging=true;lastPtr={x:e.clientX,y:e.clientY};});
    window.addEventListener('pointermove',e=>{
      if(!dragging)return;
      orbit.theta-=(e.clientX-lastPtr.x)*0.007;
      orbit.phi=Math.max(0.12,Math.min(Math.PI-0.12,orbit.phi-(e.clientY-lastPtr.y)*0.007));
      lastPtr={x:e.clientX,y:e.clientY};
    });
    window.addEventListener('pointerup',()=>{dragging=false;});
    dom.addEventListener('wheel',e=>{
      e.preventDefault();
      orbit.dist=Math.max(30,Math.min(340,orbit.dist+e.deltaY*0.12));
    },{passive:false});
  }

  function mount(container, opts){
    if(typeof THREE==='undefined') return false;
    host=container;
    if(!renderer){
      makeScene();
      renderer=new THREE.WebGLRenderer({antialias:true});
      if (ML) ML.setupRenderer(renderer);
      attach(renderer.domElement);
    }
    if(renderer.domElement.parentElement!==host) host.appendChild(renderer.domElement);
    onStats=opts.onStats||null;
    rebuild(opts);
    resize();
    running=true; raf=requestAnimationFrame(loop);
    return true;
  }

  function rebuild(opts){
    spec=buildEnzyme({
      twoDomain: !!opts.twoDomain,
      hasCleft: opts.hasCleft!==false,
      stability: opts.stability!=null?opts.stability:1,
      seed: opts.seed!=null?opts.seed:12345
    });
    eng=new Engine(spec);
    eng.kT = opts.kT!=null?opts.kT:1.2;
    buildMeshes();
    tRate=[];
    // 分子の大きさに合わせて自動でフレーミングする
    let cx=0,cy=0,cz=0;
    spec.pos.forEach(p=>{cx+=p.x;cy+=p.y;cz+=p.z;});
    cx/=spec.N; cy/=spec.N; cz/=spec.N;
    viewTarget={x:cx,y:cy,z:cz};
    let rad=0;
    spec.pos.forEach(p=>{rad=Math.max(rad,Math.hypot(p.x-cx,p.y-cy,p.z-cz));});
    orbit.dist = Math.max(55, rad*3.6);
  }

  function resize(){
    if(!renderer||!host) return;
    const w=host.clientWidth||600, h=host.clientHeight||420;
    renderer.setPixelRatio(Math.min(window.devicePixelRatio,1.6));
    renderer.setSize(w,h,false);
    camera.aspect=w/h; camera.updateProjectionMatrix();
  }

  function unmount(){ running=false; if(raf)cancelAnimationFrame(raf); raf=null; }
  function setKT(v){ if(eng) eng.kT=v; }
  function setSpeed(v){ substeps=v; }
  function reset(){ if(eng){ eng.reset(); tRate=[]; } }

  window.EnzymeMD.View = { mount, unmount, resize, setKT, setSpeed, reset, rebuild };
})();
