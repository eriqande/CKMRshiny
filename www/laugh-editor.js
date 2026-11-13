/* Discrete-age painting. SVG motion stays local; one completed stroke is sent
 * to Shiny. R owns the canonical values/history, including grid changes. */
(function () {
  "use strict";
  const NS = "http://www.w3.org/2000/svg";
  const editors = new Map();
  let latest = null, sequence = 0, registered = false;
  const bounds = {left:78, right:980, top:20, bottom:226, width:1000, height:278};
  function html(tag, attrs, text) {
    const e = document.createElement(tag);
    Object.entries(attrs || {}).forEach(([k,v]) => e.setAttribute(k,v));
    if (text !== undefined) e.textContent = text;
    return e;
  }
  function svg(tag, attrs, text) {
    const e = document.createElementNS(NS,tag);
    Object.entries(attrs || {}).forEach(([k,v])=>e.setAttribute(k,v));
    if (text !== undefined) e.textContent = text;
    return e;
  }
  function emit(name, payload) {
    if (window.Shiny) Shiny.setInputValue(name, {...payload,sequence:++sequence},{priority:"event"});
  }
  function ageLabel(index, age) { return index === age ? age+"+" : String(index); }
  function concise(value) { return Number(value.toPrecision(6)).toLocaleString(undefined,{maximumSignificantDigits:6}); }
  function niceMaximum(value) {
    value=Math.max(value,1);
    const step=Math.pow(10,Math.floor(Math.log10(value))-1);
    return Number((Math.ceil(value/step)*step).toPrecision(8));
  }
  function sendEdit(e, year=e.year, revision=latest.revision) {
    emit("laugh_curve_edit",{key:e.key, year, values:e.values.slice(), revision});
    e.status.textContent = "Edit sent. Numeric values and normalized allocations update below.";
  }
  function render(e) {
    if (!latest || !e.values) return;
    const n = e.values.length, max = Number(e.maximum.value), age = latest.plus_age;
    const step = (bounds.right-bounds.left)/n;
    const xp = i => bounds.left+step*(i+.5);
    const yp = value => bounds.bottom-Math.min(1,value/max)*(bounds.bottom-bounds.top);
    e.svg.replaceChildren();
    e.svg.appendChild(svg("text",{x:bounds.left,y:14},e.unit || ""));
    e.svg.appendChild(svg("rect",{x:bounds.left+step*age,y:bounds.top,width:step,height:bounds.bottom-bounds.top,fill:"#f8eddb"}));
    for (let t=0;t<=4;t++) {
      const value=max*t/4,y=yp(value);
      e.svg.appendChild(svg("line",{x1:bounds.left,x2:bounds.right,y1:y,y2:y,stroke:"#e2e8e4"}));
      e.svg.appendChild(svg("text",{x:bounds.left-10,y:y+4,"text-anchor":"end"},concise(value)));
    }
    const juvenile = Math.min(age-1,latest.juvenile);
    const boundary=bounds.left+step*(juvenile+1);
    e.svg.appendChild(svg("line",{x1:boundary,x2:boundary,y1:bounds.top,y2:bounds.bottom,
      stroke:"#bc423c","stroke-dasharray":"5 5","stroke-width":1.5}));
    (e.series || []).filter(s=>s.year!==e.year).forEach(s=> {
      const line=svg("polyline",{points:s.values.map((v,i)=>xp(i)+","+yp(v)).join(" "),
        fill:"none",stroke:"#aeb8ba","stroke-width":1.6,opacity:.65});
      line.appendChild(svg("title",{},s.label+" (background)")); e.svg.appendChild(line);
    });
    e.values.forEach((v,i)=> {
      e.svg.appendChild(svg("rect",{x:xp(i)-step*.33,y:yp(v),width:step*.66,height:bounds.bottom-yp(v),fill:e.color,opacity:.18}));
      if (n<=26 || i%Math.ceil(n/18)===0 || i===age)
        e.svg.appendChild(svg("text",{x:xp(i),y:bounds.bottom+22,"text-anchor":"middle"},ageLabel(i,age)));
    });
    e.svg.appendChild(svg("polyline",{points:e.values.map((v,i)=>xp(i)+","+yp(v)).join(" "),
      fill:"none",stroke:e.color,"stroke-width":e.series?3.5:2.5}));
    e.values.forEach((v,i)=>e.svg.appendChild(svg("circle",{cx:xp(i),cy:yp(v),r:3.5,fill:e.color})));
    e.svg.appendChild(svg("text",{x:(bounds.left+bounds.right)/2,y:271,"text-anchor":"middle"},"Age (final category is a plus group)"));
    e.inputs.forEach((input,i)=> { if (document.activeElement!==input) input.value=Number(e.values[i].toPrecision(10)); });
  }
  function buildTable(e) {
    const scroll=e.table.scrollLeft;
    const table=html("table"), ages=html("tr"), values=html("tr");
    ages.appendChild(html("th",{scope:"row"},"Age")); values.appendChild(html("th",{scope:"row"},"Value"));
    e.inputs=[];
    e.values.forEach((value,i)=> {
      const label=ageLabel(i,latest.plus_age);
      ages.appendChild(html("th",{scope:"col"},label));
      const cell=html("td"), input=html("input",{type:"number",min:0,step:"any",
        "aria-label":e.title+", age "+label});
      if (e.key==="pmatata") input.max=1;
      input.value=Number(value.toPrecision(10));
      input.addEventListener("change",()=> {
        const number=Number(input.value);
        if (!input.value.trim() || !Number.isFinite(number) || number<0 || (e.key==="pmatata" && number>1)) {
          input.value=Number(e.values[i].toPrecision(10));
          e.status.textContent="Enter a finite nonnegative value"+(e.key==="pmatata"?" between 0 and 1.":".");
          return;
        }
        e.values[i]=number;
        if (number>Number(e.maximum.value)) e.maximum.value=niceMaximum(number*1.1);
        render(e); sendEdit(e);
      });
      cell.appendChild(input); values.appendChild(cell); e.inputs.push(input);
    });
    table.append(ages,values); e.table.replaceChildren(table); e.table.scrollLeft=scroll;
    e.table.setAttribute("aria-label",e.title+" numeric age values; scroll sideways to older ages");
  }
  function point(e,event) {
    const local=new DOMPoint(event.clientX,event.clientY).matrixTransform(e.svg.getScreenCTM().inverse());
    const step=(bounds.right-bounds.left)/e.values.length;
    return {i:Math.max(0,Math.min(e.values.length-1,Math.floor((local.x-bounds.left)/step))),
      value:Math.max(0,Math.min(Number(e.maximum.value),
        (bounds.bottom-local.y)/(bounds.bottom-bounds.top)*Number(e.maximum.value))), x:local.x,y:local.y};
  }
  function paint(e,p) {
    const previous=e.last;
    if (previous && previous.i!==p.i) {
      const direction=p.i>previous.i?1:-1;
      for(let i=previous.i;i!==p.i;i+=direction)
        e.values[i]=previous.value+(p.value-previous.value)*(i-previous.i)/(p.i-previous.i);
    }
    e.values[p.i]=p.value; e.last=p;
    e.status.textContent="Age "+ageLabel(p.i,latest.plus_age)+": "+concise(p.value)+". Release to apply the stroke.";
    render(e);
  }
  function createEditor(root) {
    const key=root.dataset.key, e={root,key,inputs:[]};
    const controls=html("div",{class:"laugh-editor-controls"});
    const label=html("label",{},"Starting shape"), preset=html("select",{"aria-label":key+" starting shape"});
    [["example","Example"],["flat","Flat"],["increasing","Increasing"],["declining","Declining"],["bell","Bell-shaped"],["clear","Clear to zero"]]
      .forEach(([value,text])=>preset.appendChild(html("option",{value},text)));
    label.appendChild(preset); controls.appendChild(label);
    const apply=html("button",{type:"button",class:"btn btn-default"},"Apply shape");
    apply.addEventListener("click",()=>emit("laugh_curve_action",{key,year:e.year,action:"preset",preset:preset.value,revision:latest.revision}));
    controls.appendChild(apply);
    [.5,2].forEach(factor=> {
      const button=html("button",{type:"button",class:"btn btn-default","aria-label":key+" multiply all values by "+factor},factor===.5?"Values \u00d7 0.5":"Values \u00d7 2");
      button.addEventListener("click",()=>emit("laugh_curve_action",{key,year:e.year,action:"scale",factor,revision:latest.revision}));
      controls.appendChild(button);
    });
    const scaleLabel=html("label",{},"View y maximum"), maximum=html("input",{type:"number",min:.000001,step:"any","aria-label":key+" view y maximum"});
    maximum.addEventListener("change",()=> {
      if (!Number.isFinite(Number(maximum.value)) || Number(maximum.value)<=0) maximum.value=Math.max(...e.values,1);
      render(e); e.status.textContent="Display scale changed; curve values and scientific inputs are unchanged.";
    });
    scaleLabel.appendChild(maximum); controls.appendChild(scaleLabel); e.maximum=maximum;
    if (key==="pmatata") { maximum.value=1; maximum.disabled=true; }
    root.appendChild(controls);
    root.appendChild(html("p",{class:"hint"},"Drag a stroke across ages to paint values. The y maximum changes the drawing scale only. The shaded final column is the plus group."+
      (key==="pmatata"?" Maturity values, including scaled values, are capped at one.":"")));
    e.svg=svg("svg",{viewBox:"0 0 1000 278",role:"img","aria-label":key+" age curve. Draw with a pointer, or use the numeric table below."});
    root.appendChild(e.svg);
    e.status=html("p",{class:"laugh-editor-status",role:"status","aria-live":"polite"}); root.appendChild(e.status);
    const details=html("details"), summary=html("summary",{},"Edit numeric values or paste a whole curve"); details.appendChild(summary);
    details.appendChild(html("p",{class:"hint"},"Scroll sideways to see older ages and the plus group. Age and Value stay fixed on the left. You can also focus the window and use arrow keys, Home, or End."));
    e.table=html("div",{class:"laugh-editor-table",tabindex:"0",role:"region"}); details.appendChild(e.table);
    e.table.addEventListener("keydown",event=> {
      if(event.target!==e.table) return;
      const moves={ArrowLeft:-160,ArrowRight:160,Home:-e.table.scrollWidth,End:e.table.scrollWidth};
      if(Object.prototype.hasOwnProperty.call(moves,event.key)) {
        event.preventDefault(); e.table.scrollLeft+=moves[event.key];
      }
    });
    const paste=html("div",{class:"laugh-editor-paste"}), area=html("textarea",{"aria-label":key+" pasted age values",placeholder:"One value for every age, in order; final value is the whole plus group."});
    const pasteButton=html("button",{type:"button",class:"btn btn-default"},"Apply pasted values");
    pasteButton.addEventListener("click",()=> {
      const tokens=area.value.trim().split(/[,;\s]+/), numbers=tokens.map(Number);
      if (!area.value.trim() || numbers.length!==e.values.length || numbers.some(v=>!Number.isFinite(v)||v<0||(key==="pmatata"&&v>1))) {
        e.status.textContent="Paste exactly "+e.values.length+" finite nonnegative values"+(key==="pmatata"?", all at most 1":"")+" (ages 0 through "+latest.plus_age+"+)."; return;
      }
      e.values=numbers; e.maximum.value=key==="pmatata"?1:Math.max(Number(e.maximum.value),...numbers);
      render(e); sendEdit(e);
    });
    paste.append(area,pasteButton); details.appendChild(paste); root.appendChild(details);
    e.svg.addEventListener("pointerdown",event=> {
      if (event.button!==0 || !e.values) return;
      const p=point(e,event);
      if(p.x<bounds.left||p.x>bounds.right||p.y<bounds.top||p.y>bounds.bottom) return;
      event.preventDefault(); e.before=e.values.slice(); e.last=null; e.drawing=true;
      e.strokeYear=e.year; e.strokeRevision=latest.revision;
      e.svg.setPointerCapture(event.pointerId); paint(e,p);
    });
    e.svg.addEventListener("pointermove",event=> { if(e.drawing) paint(e,point(e,event)); });
    e.svg.addEventListener("pointerup",event=> {
      if(!e.drawing) return;
      paint(e,point(e,event)); e.drawing=false; e.last=null;
      if(e.svg.hasPointerCapture(event.pointerId)) e.svg.releasePointerCapture(event.pointerId);
      sendEdit(e,e.strokeYear,e.strokeRevision);
    });
    e.svg.addEventListener("pointercancel",()=> {
      if(e.drawing) { e.values=e.before; e.drawing=false; e.last=null; render(e); e.status.textContent="Stroke cancelled."; }
    });
    editors.set(key,e);
  }
  function update(message) {
    latest=message;
    if(!editors.size) return;
    message.curves.forEach(curve=> {
      const e=editors.get(curve.key); if(!e) return;
      // Never replace an in-progress stroke with a cutoff-only server update.
      if(e.drawing) return;
      const changedGrid=!e.values || e.values.length!==curve.values.length;
      const changedYear=e.year!==curve.year;
      const changedValues=changedGrid || e.values.some((value,i)=>value!==curve.values[i]);
      e.values=curve.values.slice(); e.title=curve.title; e.color=curve.color; e.unit=curve.unit;
      e.year=curve.year; e.series=curve.series;
      if(e.series) {
        e.color=e.series.find(s=>s.year===e.year).color;
        e.title+=" — year "+e.year;
        if(changedYear) e.status.textContent="Editing year "+e.year+". Other active years are muted in the background.";
      }
      const peak=Math.max(...e.values,...(e.series||[]).flatMap(s=>s.values),0);
      if(curve.key==="pmatata") e.maximum.value=1;
      else if(changedGrid || !Number(e.maximum.value) || (changedValues && peak>Number(e.maximum.value)))
        e.maximum.value=niceMaximum(peak*1.1);
      if(changedGrid || changedYear) buildTable(e);
      render(e);
    });
  }
  function ready() {
    if(window.Shiny && !registered) { Shiny.addCustomMessageHandler("laugh_editor_state",update); registered=true; }
    document.querySelectorAll(".laugh-editor").forEach(root=> {if(!editors.has(root.dataset.key)) createEditor(root);});
    if(latest) update(latest);
    emit("laugh_editor_ready",{});
  }
  if(window.Shiny) { Shiny.addCustomMessageHandler("laugh_editor_state",update); registered=true; }
  if(document.readyState==="loading") document.addEventListener("DOMContentLoaded",ready); else ready();
  if(window.jQuery) jQuery(document).on("shiny:connected",ready);
})();
