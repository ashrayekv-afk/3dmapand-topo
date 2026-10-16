/* Patient-linked corneal shape. Image estimates and measured exports are kept
 * distinct. No generic cone, fitted progression model, remote reader or PHI log. */
(function(root){'use strict';
const node=typeof module!=='undefined'&&module.exports;
const I=node?require('../image_regions.js'):root.NKPIImageRegions;
const C=node?require('./surface_schema.js'):root.CornealData;
const finite=Number.isFinite,clamp=(v,a,b)=>Math.max(a,Math.min(b,v));
function solve(a,b){a=a.map((r,i)=>[...r,b[i]]);for(let i=0;i<b.length;i++){let k=i;for(let j=i+1;j<b.length;j++)if(Math.abs(a[j][i])>Math.abs(a[k][i]))k=j;[a[i],a[k]]=[a[k],a[i]];if(Math.abs(a[i][i])<1e-10)return null;const v=a[i][i];a[i]=a[i].map(x=>x/v);for(let j=0;j<b.length;j++)if(j!==i){const f=a[j][i];a[j]=a[j].map((x,q)=>x-f*a[i][q]);}}return a.map(r=>r[b.length]);}
function smooth(grid,i,j,step){
 // Local weighted plane, including a nearest-support and coverage requirement.
 // Small printed-number occlusions may be interpolated; large gaps remain open.
 const radius=.6,n=Math.ceil(radius/step),points=[];let possible=0,near=Infinity;
 for(let y=Math.max(0,j-n);y<=Math.min(40,j+n);y++)for(let x=Math.max(0,i-n);x<=Math.min(40,i+n);x++){
  const dx=(x-i)*step,dy=(y-j)*step,d=Math.hypot(dx,dy);if(d>radius||Math.hypot(x-20,y-20)>20)continue;possible++;const q=grid[y][x];if(q){points.push({dx,dy,v:q.value,w:Math.exp(-d*d/.08)});near=Math.min(near,d);}
 }
 if(points.length<6||points.length<possible*.5||near>Math.max(step*2.01,.25))return null;
 const a=Array.from({length:3},()=>[0,0,0]),b=[0,0,0];for(const p of points){const v=[1,p.dx,p.dy];for(let k=0;k<3;k++){b[k]+=p.w*v[k]*p.v;for(let l=0;l<3;l++)a[k][l]+=p.w*v[k]*v[l];}}
 const fit=solve(a,b);return fit?{value:fit[0],dx:fit[1],dy:fit[2],interpolated:!grid[j][i]}:null;
}
function imageModel(shape,thickness){
 I.validateShape(shape);I.validateMap(thickness);
 if(thickness.id!=='thickness'||Math.abs(shape.diameter_mm-thickness.diameter_mm)>.001)throw Error('Elevation and thickness need the same printed plot diameter.');
 const radius=shape.diameter_mm/3,step=radius/20,R=shape.bfs_radius_mm,vertices=[],nx=41;let possible=0,interpolated=0;
 for(let j=0;j<41;j++)for(let i=0;i<41;i++){
  const x=(i-20)*step,y=(j-20)*step;if(Math.hypot(x,y)>radius+1e-8){vertices.push(null);continue;}possible++;
  const e=smooth(shape.elevation_um,i,j,step),t=smooth(thickness.grid,i,j,step);
  if(!e||!t||t.value<100||t.value>1200||R*R<=x*x+y*y){vertices.push(null);continue;}
  const sag=Math.sqrt(R*R-x*x-y*y),z=sag-R+e.value/1000,dx=-x/sag+e.dx/1000,dy=y/sag-e.dy/1000,den=Math.hypot(dx,dy,1),normal=[-dx/den,-dy/den,1/den];
  const front=[x,-y,z],back=front.map((v,k)=>v-normal[k]*t.value/1000);
  if(![...front,...back].every(finite)){vertices.push(null);continue;}
  if(e.interpolated||t.interpolated)interpolated++;
  vertices.push({front,back,thickness:t.value,interpolated:e.interpolated||t.interpolated});
 }
 const valid=vertices.filter(Boolean);if(valid.length/possible<.5)throw Error('Too little overlapping elevation and thickness coverage for a 3D estimate.');
 return finish({vertices,nx,ny:41,kind:'estimated',radius,coverage:valid.length/possible,interpolated,possible,bfs:R,fit_diameter:shape.fit_diameter_mm,legend_step:shape.legend_step_um});
}
function nativeModel(raw,eye,visit){
 const s=C.validate(raw);if(s.kind!=='measured'||s.eye!==eye||s.visit!==visit||!s.anterior?.z_mm)throw Error('A matching measured anterior height grid is needed.');
 const cs=s.coordinate_system,mx=(cs.x_positive==='temporal'?1:-1)*(eye==='OD'?-1:1),my=cs.y_positive==='superior'?1:-1,mz=cs.z_positive==='anterior'?1:-1;
 const vertices=[];let possible=0;for(let j=0;j<s.y_mm.length;j++)for(let i=0;i<s.x_mm.length;i++){
  const x=s.x_mm[i]*mx,y=s.y_mm[j]*my;if(Math.hypot(x,y)>3){vertices.push(null);continue;}possible++;
  const z=s.anterior.z_mm[j][i],b=s.posterior?.z_mm?.[j]?.[i],t=s.thickness_um?.[j]?.[i];
  if(!finite(z)||(finite(b)&&mz*z<=mz*b))vertices.push(null);else vertices.push({front:[x,y,mz*z],back:finite(b)?[x,y,mz*b]:null,thickness:finite(t)?t:null,interpolated:false});
 }
 if(!vertices.some(Boolean))throw Error('No usable height samples in the central region.');
 return finish({vertices,nx:s.x_mm.length,ny:s.y_mm.length,kind:'native',radius:3,coverage:vertices.filter(Boolean).length/possible,possible,interpolated:0});
}
function finish(m){
 const triangles=[];for(let j=0;j<m.ny-1;j++)for(let i=0;i<m.nx-1;i++){const a=j*m.nx+i,b=a+1,c=a+m.nx,d=c+1;for(const ids of [[a,b,c],[b,d,c]])if(ids.every(k=>m.vertices[k]))triangles.push(ids);}
 if(!triangles.length)throw Error('No complete surface cells.');m.triangles=triangles;
 const ts=m.vertices.flatMap(v=>v&&finite(v.thickness)?[v.thickness]:[]);m.thicknessRange=ts.length?[Math.min(...ts),Math.max(...ts)]:null;
 m.centerZ=m.vertices.filter(Boolean).reduce((s,v)=>s+v.front[2],0)/m.vertices.filter(Boolean).length;
 return m;
}
function fromReport(s){
 if(s.surface?.kind==='measured'&&s.surface?.anterior?.z_mm)return nativeModel(s.surface,s.eye,s.visit);
 if(!s.shape_estimate)throw Error(s.shape_issue||'A readable Belin–Ambrósio elevation and thickness page is needed in this selected study.');
 const ts=(s.image_maps||[]).filter(v=>v.id==='thickness');if(ts.length!==1)throw Error('One unambiguous thickness map is needed.');
 return imageModel(s.shape_estimate,ts[0]);
}
const api={smooth,imageModel,nativeModel,fromReport};if(node){module.exports=api;return;}
const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function mount(el,sets,eye){
 let visit=sets.find(s=>s.visit==='latest')?'latest':'first',yaw=.25,pitch=1.08,cut=false,color='thickness',model=null,error='',drag=null,dead=false;
 const data=()=>sets.find(s=>s.visit===visit);
 function load(){try{model=fromReport(data());error='';}catch(e){model=null;error=e.message;}}
 function layout(){
 load();const s=data();
 el.innerHTML='<div class="rg-shape-toolbar"><div class="rg-actions">'+sets.map(v=>'<button type="button" data-shape-visit="'+v.visit+'" aria-pressed="'+(visit===v.visit)+'">'+(v.visit==='latest'?'Latest scan':sets.length>1?'Earlier scan':'Selected scan')+'</button>').join('')+'</div>'+(model?'<div class="rg-actions"><button type="button" data-shape-view="front">Front</button><button type="button" data-shape-view="profile">Profile</button><button type="button" data-shape-view="reset">Reset view</button></div>':'')+'</div>';
 if(!model){el.innerHTML+='<p class="rg-empty-small">'+esc(error)+' The original report remains available below.</p>';return;}
 const estimate=model.kind==='estimated';
 el.innerHTML+='<div class="rg-shape-label"><b>'+eye+' · '+(estimate?'Image-derived corneal estimate':'Device height reconstruction')+'</b><span>'+Math.round(model.coverage*100)+'% modeled coverage · central '+(2*model.radius).toFixed(1)+' mm</span></div><p class="rg-small">'+(estimate?'Anterior shape from elevation + printed BFS. Inner surface approximated from thickness; the floating reference’s position and tilt are not available.':'Front and back use supplied height coordinates. Thickness colors use device pachymetry, when supplied.')+'</p><div class="rg-shape-stage"><canvas width="1000" height="650" tabindex="0" aria-label="'+eye+' 3D corneal reconstruction. Drag or use arrow keys to rotate; Front and Profile buttons reset orientation."></canvas></div><div class="rg-shape-toolbar"><label>Color <select data-shape-color><option value="thickness">Thickness</option><option value="shape">Surface shape</option></select></label><label><input type="checkbox" data-shape-cut> Cutaway</label><button type="button" data-shape-save>Save 3D figure</button></div><div class="rg-shape-key">'+(model.thicknessRange?'<span>Thinner</span><i></i><span>Thicker</span><b>~'+model.thicknessRange.map(v=>Math.round(v)).join('–')+' µm</b>':'No thickness overlay supplied.')+'</div><p class="rg-small">Drag to rotate · same scale in all directions · open patches mean missing image data, not tissue defects. '+(estimate?'Small text occlusions are locally interpolated. This is an approximate corneal segment, not an exact replica or a clinical measurement surface.':'The corneal segment is limited to the supplied grid; it is not a full-eye model.')+'</p><details class="rg-details"><summary>Reconstruction details</summary><p>'+ (estimate?'BFS radius '+model.bfs.toFixed(2)+' mm; fitting diameter '+(model.fit_diameter===null?'unreadable':model.fit_diameter.toFixed(2))+' mm. The centered sphere’s sag is combined with signed elevation (µm → mm). A 0.6-mm local weighted plane smooths readable color samples; only small gaps with nearby support are interpolated. The back surface follows the anterior normal using estimated pachymetry. This does not recover the separately measured posterior surface.':'Only native z coordinates determine the geometry. Native grid orientation is converted for this eye. Surfaces with reversed front/back ordering are omitted.')+'</p><p>Visits are reconstructed independently; switching them does not establish a registered height change. Use the regional maps for the image comparison. NKPI is unchanged.</p><p>'+model.interpolated+' of '+model.possible+' sample positions required local interpolation. Review the original map if quality is poor or the reconstruction is unexpected.</p></details>';
 el.querySelector('[data-shape-cut]').checked=cut;el.querySelector('[data-shape-color]').value=color;
 const canvas=el.querySelector('canvas');canvas.addEventListener('pointerdown',e=>{drag={x:e.clientX,y:e.clientY};canvas.setPointerCapture(e.pointerId);});canvas.addEventListener('pointermove',e=>{if(!drag)return;yaw+=(e.clientX-drag.x)*.008;pitch=clamp(pitch+(e.clientY-drag.y)*.008,-1.5,1.5);drag={x:e.clientX,y:e.clientY};draw();});canvas.addEventListener('pointerup',()=>drag=null);canvas.addEventListener('pointercancel',()=>drag=null);canvas.addEventListener('keydown',e=>{if(!['ArrowLeft','ArrowRight','ArrowUp','ArrowDown'].includes(e.key))return;e.preventDefault();if(e.key==='ArrowLeft')yaw-=.1;if(e.key==='ArrowRight')yaw+=.1;if(e.key==='ArrowUp')pitch=clamp(pitch-.1,-1.5,1.5);if(e.key==='ArrowDown')pitch=clamp(pitch+.1,-1.5,1.5);draw();});draw();
 }
 function project(v){const z=v[2]-model.centerZ,x=v[0]*Math.cos(yaw)+z*Math.sin(yaw),z1=-v[0]*Math.sin(yaw)+z*Math.cos(yaw),y=v[1]*Math.cos(pitch)-z1*Math.sin(pitch),depth=v[1]*Math.sin(pitch)+z1*Math.cos(pitch);const k=300/model.radius;return [500+x*k,315-y*k,depth];}
 function tint(t,back){if(back)return '#87b4b0';if(color==='shape'||!model.thicknessRange||!finite(t))return '#a18ac7';const [lo,hi]=model.thicknessRange,f=clamp((t-lo)/Math.max(1,hi-lo),0,1),a=[115,69,158],b=[204,230,226];return 'rgb('+a.map((v,k)=>Math.round(v+(b[k]-v)*f)).join(',')+')';}
 function draw(){if(dead||!model)return;const c=el.querySelector('canvas'),ctx=c.getContext('2d');ctx.clearRect(0,0,1000,650);ctx.fillStyle='#faf8fc';ctx.fillRect(0,0,1000,650);const faces=[];
  for(const ids of model.triangles){const vs=ids.map(i=>model.vertices[i]);if(cut&&vs.some(v=>v.front[0]>1e-8))continue;for(const side of ['back','front']){if(vs.some(v=>!v[side]))continue;const points=vs.map(v=>project(v[side]));faces.push({points,fill:tint(vs.reduce((s,v)=>s+(v.thickness||0),0)/3,side==='back'),depth:points.reduce((s,v)=>s+v[2],0)/3});}
  }
  // Boundary walls join only existing front/back samples, including cut edges.
  const edges=new Map();for(const ids of model.triangles){if(cut&&ids.some(i=>model.vertices[i].front[0]>1e-8))continue;for(const [a,b]of [[ids[0],ids[1]],[ids[1],ids[2]],[ids[2],ids[0]]]){const key=[a,b].sort((a,b)=>a-b).join(':');if(edges.has(key))edges.delete(key);else edges.set(key,[a,b]);}}
  for(const [a,b]of edges.values()){const v=model.vertices[a],w=model.vertices[b];if(!v.back||!w.back)continue;const points=[v.front,w.front,w.back,v.back].map(project);faces.push({points,fill:'#cbbdde',depth:points.reduce((s,p)=>s+p[2],0)/4});}
  faces.sort((a,b)=>a.depth-b.depth);for(const f of faces){ctx.beginPath();f.points.forEach((p,i)=>i?ctx.lineTo(p[0],p[1]):ctx.moveTo(p[0],p[1]));ctx.closePath();ctx.fillStyle=f.fill;ctx.fill();ctx.strokeStyle=f.fill;ctx.lineWidth=.6;ctx.stroke();}
  ctx.font='18px system-ui';ctx.textAlign='center';ctx.fillStyle='#514461';const r=model.radius*1.12;for(const [label,v]of [['Superior',[0,r,0]],['Inferior',[0,-r,0]],[eye==='OD'?'Temporal':'Nasal',[-r,0,0]],[eye==='OD'?'Nasal':'Temporal',[r,0,0]]]){const p=project(v);ctx.fillText(label,p[0],p[1]);}
  ctx.textAlign='left';ctx.font='bold 20px system-ui';ctx.fillText(eye+' · '+(visit==='latest'?'Latest':'Earlier / selected')+' · '+(model.kind==='estimated'?'IMAGE-DERIVED ESTIMATE':'DEVICE HEIGHTS'),24,30);ctx.font='16px system-ui';ctx.fillText('True proportions · '+(model.kind==='estimated'?'inner surface inferred from thickness':'native coordinates'),24,625);
 }
 el.addEventListener('click',click);el.addEventListener('change',change);
 function click(e){const b=e.target.closest('button');if(!b)return;if(b.dataset.shapeVisit){visit=b.dataset.shapeVisit;layout();}else if(b.dataset.shapeView){yaw=b.dataset.shapeView==='reset'?.25:0;pitch=b.dataset.shapeView==='front'?0:b.dataset.shapeView==='profile'?1.48:1.08;draw();}else if(b.hasAttribute('data-shape-save')&&model){const canvas=el.querySelector('canvas'),a=document.createElement('a');a.href=canvas.toDataURL('image/png');a.download='NKPI_'+eye+'_'+visit+'_3D_'+model.kind+'.png';a.click();}}
 function change(e){if(e.target.hasAttribute('data-shape-cut'))cut=e.target.checked;else if(e.target.hasAttribute('data-shape-color'))color=e.target.value;else return;draw();}
 layout();return ()=>{dead=true;el.removeEventListener('click',click);el.removeEventListener('change',change);};
}
root.NKPIPatientShape={...api,mount};
})(typeof window!=='undefined'?window:globalThis);
