/*
 * 分子アニメーション共通の見た目（Three.js）
 *
 * 分子生物学の3Dアニメーション（WEHI の Drew Berry らの作品など）に見られる
 * 「暗い背景・霧による奥行き・輪郭が光るリムライト・でこぼこした分子表面」を
 * 3つのビュー（転写・翻訳 / 酵素 / 細胞内夾雑）で共通に使うための小道具です。
 */
(function(){
  'use strict';
  if (typeof THREE === 'undefined') return;

  const BG = 0x050b14;

  // sRGB で指定した色を、ライティング計算用の線形色に直す
  function lin(hex){ return new THREE.Color(hex).convertSRGBToLinear(); }

  // MeshStandardMaterial に「視線と面が浅い角度になるほど光る」リム項を足す。
  // 定数はシェーダに埋め込まず uniform で渡す（three.js はシェーダを
  // onBeforeCompile の関数本体でキャッシュするため、埋め込むと別マテリアルと混ざる）。
  function material(hex, o){
    o = o || {};
    const mat = new THREE.MeshStandardMaterial({
      color: o.vertexColor ? 0xffffff : lin(hex),
      roughness: o.roughness != null ? o.roughness : 0.62,
      metalness: 0.0,
      emissive: o.emissive != null ? lin(o.emissive) : new THREE.Color(0x000000),
      emissiveIntensity: o.emissiveIntensity != null ? o.emissiveIntensity : 1,
      transparent: !!o.transparent,
      opacity: o.opacity != null ? o.opacity : 1,
      depthWrite: o.depthWrite != null ? o.depthWrite : true
    });
    const rim = { value: lin(o.rim != null ? o.rim : 0x9ad7ff) };
    const rimPow = { value: o.rimPow != null ? o.rimPow : 2.4 };
    const rimStr = { value: o.rimStr != null ? o.rimStr : 0.55 };
    mat.userData.rimStr = rimStr;
    mat.onBeforeCompile = sh => {
      sh.uniforms.uRim = rim; sh.uniforms.uRimPow = rimPow; sh.uniforms.uRimStr = rimStr;
      sh.fragmentShader = 'uniform vec3 uRim;\nuniform float uRimPow;\nuniform float uRimStr;\n' +
        sh.fragmentShader.replace('#include <tonemapping_fragment>',
          'float rimF = pow(1.0 - clamp(abs(dot(normalize(normal), normalize(vViewPosition))), 0.0, 1.0), uRimPow);\n' +
          'gl_FragColor.rgb += uRim * rimF * uRimStr;\n' +
          '#include <tonemapping_fragment>');
    };
    return mat;
  }

  // レンダラーの色・トーンの設定
  function setupRenderer(renderer){
    renderer.outputEncoding = THREE.sRGBEncoding;
    renderer.toneMapping = THREE.ACESFilmicToneMapping;
    renderer.toneMappingExposure = 1.0;
  }

  // 暗い背景＋霧＋3点照明（上からの柔らかい光、キーライト、背後からの青いリム光）
  function setupScene(scene, fogDensity){
    scene.background = new THREE.Color(BG);
    scene.fog = new THREE.FogExp2(BG, fogDensity != null ? fogDensity : 0.0032);
    scene.add(new THREE.HemisphereLight(0xcfe6ff, 0x1a1020, 0.75));
    const key = new THREE.DirectionalLight(0xfff1dc, 1.25); key.position.set(0.6, 1.0, 0.8); scene.add(key);
    const fill = new THREE.DirectionalLight(0x7aa7ff, 0.45); fill.position.set(-1, -0.4, 0.3); scene.add(fill);
    const back = new THREE.DirectionalLight(0x66c8ff, 0.9); back.position.set(-0.3, 0.4, -1); scene.add(back);
    return { key, fill, back };
  }

  // 遠景にぼんやり浮かぶ「ほかの分子たち」。霧に溶けて奥行きを出す。
  function Backdrop(scene, count, spread){
    this.spread = spread;
    this.items = [];
    const palette = [0x2d5a6e, 0x4a3b6b, 0x2f5a46, 0x5a4636, 0x3a4f7a];
    const geo = new THREE.IcosahedronGeometry(1, 2);
    const mat = material(0xffffff, { vertexColor: true, rimStr: 0.22, roughness: 0.85 });
    this.mesh = new THREE.InstancedMesh(geo, mat, count);
    this.mesh.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(count * 3), 3);
    this.mesh.frustumCulled = false;
    this.mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    const c = new THREE.Color();
    for (let i = 0; i < count; i++) {
      const it = {
        off: new THREE.Vector3((Math.random()-0.5)*2, (Math.random()-0.5)*2, (Math.random()-0.5)*2).multiplyScalar(spread),
        r: 8 + Math.random() * 18,
        v: new THREE.Vector3()
      };
      if (it.off.length() < spread * 0.45) it.off.setLength(spread * (0.45 + Math.random() * 0.5));
      this.items.push(it);
      c.copy(lin(palette[i % palette.length])); this.mesh.setColorAt(i, c);
    }
    scene.add(this.mesh);
    this._m = new THREE.Matrix4(); this._q = new THREE.Quaternion(); this._s = new THREE.Vector3();
    this._p = new THREE.Vector3();
  }
  // 分子は空間に固定してゆっくり漂わせ、視点(center)から離れすぎたものは
  // 反対側へ回り込ませる（カメラが動くと遠景がゆっくり流れて見える）
  // camera を渡すと、注目点(center)より手前にある分子は描かない（視界をふさがないように）
  Backdrop.prototype.update = function(center, dt, camera){
    const S = this.spread;
    const camD = camera ? camera.position.distanceTo(center) : 0;
    this.items.forEach((it, i) => {
      it.v.x += (Math.random()-0.5) * 6 * dt; it.v.y += (Math.random()-0.5) * 6 * dt; it.v.z += (Math.random()-0.5) * 6 * dt;
      it.v.multiplyScalar(0.98);
      it.off.addScaledVector(it.v, dt * 4);
      ['x','y','z'].forEach(a => {
        const d = it.off[a] - center[a];
        if (d > S) it.off[a] -= 2*S;
        if (d < -S) it.off[a] += 2*S;
      });
      this._p.copy(it.off);
      const behind = !camera || camera.position.distanceTo(this._p) > camD + it.r + 40;
      const r = behind ? it.r : 0;
      this._s.set(r, r, r);
      this._m.compose(this._p, this._q, this._s);
      this.mesh.setMatrixAt(i, this._m);
    });
    this.mesh.instanceMatrix.needsUpdate = true;
  };

  window.MolLook = { BG, lin, material, setupRenderer, setupScene, Backdrop };
})();
