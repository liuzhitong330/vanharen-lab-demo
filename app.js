(function(){
  var D = window.VH;
  var LABELS = D.labels, PCD = D.pcd, BIRTH = D.birth_pcd;
  var GENES = D.genes, RAW = D.raw, HL = D.hl;
  var NAGE = LABELS.length;
  var SVGNS = "http://www.w3.org/2000/svg";

  var CLASS = {
    myelin:{name:"a myelin gene", color:"#2f7d5f"},
    perox:{name:"a peroxisomal / very-long-chain-fat gene", color:"#b0563e"},
    lipid:{name:"a lipid or cholesterol enzyme", color:"#8a6d3b"},
    immune:{name:"a microglia / immune gene", color:"#2f6f9f"},
    other:{name:"", color:"#5a6b7a"}
  };

  // childhood window = birth (280 pcd) to 8 years (280 + 8*365.25)
  var CH_LO = BIRTH, CH_HI = BIRTH + 8*365.25;
  var chIdx=[], preIdx=[];
  for(var i=0;i<NAGE;i++){
    if(PCD[i]>=CH_LO && PCD[i]<=CH_HI) chIdx.push(i);
    if(PCD[i]<BIRTH) preIdx.push(i);
  }
  // birth boundary between last prenatal and first postnatal
  var birthAt = preIdx.length; // fractional index position of birth (between preIdx.last and next)

  function el(tag, attrs, txt){
    var e=document.createElementNS(SVGNS, tag);
    for(var k in attrs) e.setAttribute(k, attrs[k]);
    if(txt!=null) e.textContent=txt;
    return e;
  }
  function colorFor(g){ var c=HL[g]||"other"; return (CLASS[c]||CLASS.other).color; }
  function esc(s){return (s+"").replace(/[&<>]/g,function(c){return{"&":"&amp;","<":"&lt;",">":"&gt;"}[c];});}

  // ---- autocomplete ----
  var names = Object.keys(GENES).sort();
  (function(){
    var dl=document.getElementById("genelist"), frag=document.createDocumentFragment();
    for(var i=0;i<names.length;i++){var o=document.createElement("option");o.value=names[i];frag.appendChild(o);}
    dl.appendChild(frag);
  })();
  function findCI(name){var up=name.toUpperCase();for(var i=0;i<names.length;i++){if(names[i].toUpperCase()===up)return names[i];}return null;}

  // ---- trajectory plot (viewBox 0 0 600 380) ----
  var W=600, H=380, L=52, R=18, T=44, B=60;
  var iw=W-L-R, ih=H-T-B;
  function X(i){ return L + (NAGE===1?iw/2 : i*iw/(NAGE-1)); }

  // tick labels to show (avoid crowding)
  var SHOWTICK = {};
  ["8pcw","19pcw","4mo","1y","4y","8y","19y","40y"].forEach(function(l){ SHOWTICK[l]=1; });

  function draw(gene){
    var m=GENES[gene].m, raw=RAW[gene]||null, col=colorFor(gene);
    var svg=document.getElementById("scatter"); svg.innerHTML="";

    var lo=99, hi=-99;
    for(var i=0;i<m.length;i++){ if(m[i]<lo)lo=m[i]; if(m[i]>hi)hi=m[i]; }
    if(raw) raw.forEach(function(a){a.forEach(function(v){if(v<lo)lo=v;if(v>hi)hi=v;});});
    if(hi-lo<1){hi+=0.5;lo-=0.5;}
    var pad=(hi-lo)*0.12; lo=Math.max(0,lo-pad); hi+=pad;
    function Y(v){return T + ih*(1-(v-lo)/(hi-lo));}

    // childhood window shading
    if(chIdx.length){
      var x0=X(chIdx[0]), x1=X(chIdx[chIdx.length-1]);
      svg.appendChild(el("rect",{x:x0,y:T,width:(x1-x0),height:ih,fill:"#2f7d5f",opacity:0.07}));
      svg.appendChild(el("text",{x:(x0+x1)/2,y:T-24,"text-anchor":"middle","font-size":11,"font-family":"system-ui,sans-serif",fill:"#2f7d5f","font-weight":600},"early childhood"));
      svg.appendChild(el("text",{x:(x0+x1)/2,y:T-11,"text-anchor":"middle","font-size":10,"font-family":"system-ui,sans-serif",fill:"#2f7d5f"},"ALD risk window"));
    }
    // birth line (between last prenatal and first postnatal index)
    var bx = (X(birthAt-1)+X(birthAt))/2;
    svg.appendChild(el("line",{x1:bx,y1:T-2,x2:bx,y2:T+ih,stroke:"#999","stroke-dasharray":"3 3"}));
    svg.appendChild(el("text",{x:bx,y:T+ih+34,"text-anchor":"middle","font-size":9.5,"font-family":"system-ui,sans-serif",fill:"#999"},"birth"));

    // y grid
    var yt=niceTicks(lo,hi,4);
    yt.forEach(function(t){
      svg.appendChild(el("line",{x1:L,y1:Y(t),x2:L+iw,y2:Y(t),stroke:"#eee"}));
      svg.appendChild(el("text",{x:L-8,y:Y(t)+3.5,"text-anchor":"end","font-size":10.5,"font-family":"system-ui,sans-serif",fill:"#999"},t.toFixed(0)));
    });
    svg.appendChild(el("text",{x:14,y:T+ih/2,"font-size":11,"font-family":"system-ui,sans-serif",fill:"#888","text-anchor":"middle",transform:"rotate(-90 14 "+(T+ih/2)+")"},"log₂ RPKM"));

    // individual dots
    if(raw){
      for(var i=0;i<raw.length;i++){
        var arr=raw[i], cx=X(i);
        for(var j=0;j<arr.length;j++){
          var jx=cx + ((j%2?1:-1)*(1.5+(j*1.7)%5));
          svg.appendChild(el("circle",{cx:jx,cy:Y(arr[j]),r:2,fill:col,opacity:0.32}));
        }
      }
    }
    // mean line
    var d="";
    for(var i=0;i<m.length;i++){ d+=(i?"L":"M")+X(i).toFixed(1)+" "+Y(m[i]).toFixed(1)+" "; }
    svg.appendChild(el("path",{d:d,fill:"none",stroke:col,"stroke-width":2.4}));
    // markers
    for(var i=0;i<m.length;i++){ svg.appendChild(el("circle",{cx:X(i),cy:Y(m[i]),r:2.6,fill:col})); }

    // x tick labels
    for(var i=0;i<NAGE;i++){
      if(SHOWTICK[LABELS[i]]){
        svg.appendChild(el("line",{x1:X(i),y1:T+ih,x2:X(i),y2:T+ih+4,stroke:"#bbb"}));
        svg.appendChild(el("text",{x:X(i),y:T+ih+18,"text-anchor":"middle","font-size":10,"font-family":"system-ui,sans-serif",fill:"#777"},LABELS[i]));
      }
    }
  }

  function niceTicks(lo,hi,n){
    var span=hi-lo, step=Math.pow(10,Math.floor(Math.log10(span/n)));
    var err=(span/n)/step; if(err>=5)step*=5;else if(err>=2)step*=2;
    var out=[],start=Math.ceil(lo/step)*step;
    for(var v=start;v<=hi+1e-9;v+=step)out.push(Math.round(v/step)*step);
    return out;
  }

  function renderLegend(gene){
    var col=colorFor(gene), raw=RAW[gene];
    document.getElementById("legend").innerHTML=
      '<span style="--c:'+col+'">expression across development</span>'
      + (raw?'<span style="--c:'+col+';opacity:.4">individual samples</span>':'')
      + '<span style="--c:#2f7d5f;opacity:.25">early-childhood window</span>';
  }

  function describe(gene){
    var m=GENES[gene].m;
    var peakIdx=0,peakV=m[0];
    for(var i=1;i<m.length;i++){if(m[i]>peakV){peakV=m[i];peakIdx=i;}}
    // fetal baseline vs childhood peak
    var pre=0; preIdx.forEach(function(i){pre+=m[i];}); pre/=preIdx.length;
    var chPeak=-99; chIdx.forEach(function(i){if(m[i]>chPeak)chPeak=m[i];});
    var rise=chPeak-pre;
    var s="<b>"+gene+"</b> is highest at <b>"+LABELS[peakIdx]+"</b>. ";
    var inWindow = (PCD[peakIdx]>=CH_LO && PCD[peakIdx]<=CH_HI);
    if(rise>0.6){
      s+="It climbs about "+rise.toFixed(1)+" log&#8322; units from the fetal baseline into the early-childhood window. ";
    } else if(rise< -0.6){
      s+="It is actually lower in early childhood than before birth. ";
    } else {
      s+="It stays fairly flat across the childhood window. ";
    }
    var cl=HL[gene];
    if(cl==="myelin"){ s+="As "+CLASS.myelin.name+", its rise traces the myelination the young brain pays for. ";}
    else if(cl==="perox"){ s+="It is "+CLASS.perox.name+", the pathway that fails in ALD. ";}
    else if(cl==="lipid"){ s+="It is "+CLASS.lipid.name+" feeding myelin membrane. ";}
    else if(cl==="immune"){ s+="It is "+CLASS.immune.name+". ";}
    if(inWindow){ s+="Its peak lands inside the window when ALD tends to strike."; }
    else { s+="Compare its shape with MBP to see how it sits against the myelination surge."; }
    return s;
  }

  function updateMetric(gene){
    var m=GENES[gene].m, peakIdx=0,peakV=m[0];
    for(var i=1;i<m.length;i++){if(m[i]>peakV){peakV=m[i];peakIdx=i;}}
    var inWindow=(PCD[peakIdx]>=CH_LO && PCD[peakIdx]<=CH_HI);
    document.getElementById("m-dir").textContent = inWindow ? "early childhood"
        : (PCD[peakIdx]<BIRTH ? "before birth" : "adulthood");
    document.getElementById("m-dir-d").textContent = gene+" peaks here";
  }

  function pick(nameRaw){
    var name=(nameRaw||"").trim(); if(!name) return;
    var key = GENES[name] ? name : findCI(name);
    var read=document.getElementById("readout");
    if(!key){
      document.getElementById("scatter").innerHTML="";
      document.getElementById("legend").innerHTML="";
      read.innerHTML='<b>'+esc(name)+'</b> is not in this explorer. Try an official gene symbol, for example MBP, ABCD1, ACOX1, HMGCR, PLP1 or TREM2.';
      return;
    }
    document.getElementById("search").value=key;
    draw(key); renderLegend(key); updateMetric(key);
    read.innerHTML=describe(key);
  }
  window.pick=pick;

  // ---- Analysis 2: early-childhood surge ranking ----
  (function(){
    // drop non-coding / pseudogene / clone-named entries so the biology reads clearly
    var JUNK=/(\.|-|^RP\d|^RN\d|^RNU|^RNA\d|^LINC|^MIR|^SNOR|^MT-|^MTND|^MTRNR|^MTCO|^MTATP|^Z\d|^XX|^BX\d)/;
    var arr=[];
    for(var g in GENES){
      if(JUNK.test(g)) continue;
      var m=GENES[g].m;
      var pre=0; preIdx.forEach(function(i){pre+=m[i];}); pre/=preIdx.length;
      var chPeak=-99; chIdx.forEach(function(i){if(m[i]>chPeak)chPeak=m[i];});
      if(chPeak<3) continue;           // require a real childhood level
      arr.push([g, chPeak-pre]);
    }
    arr.sort(function(a,b){return b[1]-a[1];});
    var top=arr.slice(0,15), mx=top[0][1];
    var box=document.getElementById("bars");
    top.forEach(function(t){
      var cl=HL[t[0]];
      var col = cl==="myelin"?"#2f7d5f":(cl==="perox"?"#b0563e":(cl==="lipid"?"#8a6d3b":(cl==="immune"?"#2f6f9f":"#7a8794")));
      var row=document.createElement("div"); row.className="brow";
      row.innerHTML='<div class="bname">'+t[0]+'</div>'
        +'<div class="btrack"><div class="bfill" style="width:'+(100*t[1]/mx).toFixed(1)+'%;background:'+col+'"></div></div>'
        +'<div class="bnum">+'+t[1].toFixed(1)+' log₂</div>';
      row.onclick=(function(n){return function(){pick(n);window.scrollTo({top:0,behavior:"smooth"});};})(t[0]);
      box.appendChild(row);
    });
  })();

  // ---- search wiring ----
  var input=document.getElementById("search");
  input.addEventListener("change",function(){pick(input.value);});
  input.addEventListener("keydown",function(e){if(e.key==="Enter")pick(input.value);});

  pick("MBP");
})();
