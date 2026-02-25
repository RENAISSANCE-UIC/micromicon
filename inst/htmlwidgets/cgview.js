HTMLWidgets.widget({

  name: 'cgview',
  type: 'output',

  factory: function(el, width, height) {

    var cgv = null;

    return {

      renderValue: function(x) {
        // Tear down any previous instance
        el.innerHTML = '';
        cgv = null;

        // CGView.Viewer requires a CSS selector — ensure the element has an id
        if (!el.id) {
          el.id = 'cgview_' + Math.random().toString(36).substring(2, 9);
        }

        var w = el.offsetWidth  || width  || 2000;
        var h = el.offsetHeight || height || 2000;

        try {
          cgv = _render(el, w, h, x.cgview_json);
        } catch(err) {
          el.innerHTML =
            '<div style="font:12px/1.6 monospace;color:#c0392b;padding:14px;">' +
            '<strong>CGView init error:</strong><br>' + err.message + '</div>';
          console.error('[beta_cgview] init error:', err);
        }
      },

      resize: function(w2, h2) {
        if (!cgv) return;
        try { cgv.canvas.resize(w2, h2); } catch(e) {
          try { cgv.draw(); } catch(e2) { /* silent */ }
        }
      }

    };

    // =========================================================================
    // Main renderer
    // =========================================================================

    function _render(el, w, h, payload) {
      // -----------------------------------------------------------------
      // Plane C — Viewer constructor: dimensions only.
      // Everything else (backbone, ruler, feature layout) lives either in
      // the JSON (Plane A/B) or is set via JS after loadJSON (Plane C cont).
      // -----------------------------------------------------------------
      var viewer = new CGView.Viewer('#' + el.id, {
        width:  w,
        height: h
      });

      // -----------------------------------------------------------------
      // Plane A/B — Load the CGView-native JSON produced by R.
      // The payload contains: settings, sequence (length only → no GC
      // auto-tracks), backbone, ruler, features, legend, tracks.
      // -----------------------------------------------------------------
      viewer.io.loadJSON(payload);

      // -----------------------------------------------------------------
      // Plane C (continued) — post-loadJSON JS overrides, pre-draw.
      // -----------------------------------------------------------------

      // 1. arrowHeadLength: CGView's default is near-invisible (≈0.05).
      //    Set it here as a JS backstop even though it is also in the JSON
      //    settings, because loadJSON normalisation sometimes clobbers it.
      try {
        if (viewer.settings) viewer.settings.arrowHeadLength = 0.3;
      } catch(e) {}

      // 1b. Inner ruler fix.
      //   ruler.draw(innerR, outerR) immediately does innerR -= this.spacing
      //   before drawing, which pushes the inner ruler ring to a negative
      //   (off-canvas) radius when spacing is large (≥ ~500 kbp arc equivalent).
      //   Setting rulerPadding = 0 minimises the pre-draw inward offset so
      //   the inner ring lands near the edge of the innermost track rather
      //   than behind the canvas.  A longer tickLength makes the marks visible.
      try {
        if (viewer.ruler) {
          // spacing is a PIXEL gap (default 2px) — never set it to a bp value
          // or the ruler rings fly millions of pixels off-canvas.
          // tickCount (in the JSON) drives d3.ticks() → nice Mbp/kbp labels.
          // font is set correctly by loadJSON — direct string assignment here
          // would replace the Font object and break this.font.css in draw().
          // tickFormater is a getter-only computed in _updateTicks(); CGView's
          // drawLabel already appends " bp" to d3 prefix output ("1M" → "1 Mbp").
          viewer.ruler.rulerPadding = 1;
          viewer.ruler.tickLength   = 10;
          viewer.ruler.tickWidth    = 2;

          console.log('[beta_cgview] ruler after load:',
            'visible=',    viewer.ruler.visible,
            'spacing=',    viewer.ruler.spacing,
            'tickCount=',  viewer.ruler.tickCount,
            'tickLength=', viewer.ruler.tickLength,
            'font=',       viewer.ruler.font
          );
        }
      } catch(e) {
        console.warn('[beta_cgview] ruler adjustment failed:', e);
      }

      // 2. Arrow decoration.
      //    In CGView 1.8, decoration is driven by the *legend item*, not the
      //    feature itself (same pattern as swatchColor drives fill colour).
      //    The JSON already sets decoration:"arrow" on the CDS legend item;
      //    this JS pass is a backstop in case loadJSON normalisation drops it.
      var arrowLegends = { 'CDS': true, 'gene': true };
      var cgvLegend = viewer.legend;
      if (cgvLegend) {
        var legendItems = (typeof cgvLegend.items === 'function')
                            ? cgvLegend.items()
                            : (cgvLegend.items || []);
        legendItems.forEach(function(item) {
          if (arrowLegends[item.name]) {
            try { item.decoration = 'arrow'; } catch(e) {}
          }
        });
      }

      // 3. Per-track thicknessRatio.
      //    Must also be after loadJSON, before draw().
      var ratios = {
        'CDS':      2.0,   // dominant feature track — thick arcs + arrows
        'ncRNA':    0.75,  // tRNA / rRNA / repeat_region — thinner band
        'Variants': 2.0    // variant markers — same weight as CDS
      };

      var cgvTracks = (typeof viewer.tracks === 'function')
                        ? viewer.tracks()
                        : (viewer.tracks || []);

      cgvTracks.forEach(function(t) {
        var r = ratios[t.name];
        if (typeof r === 'number') {
          try { t.thicknessRatio = r; } catch(e) {
            console.warn('[beta_cgview] could not set thicknessRatio on', t.name, e);
          }
        }
      });

      viewer.draw();

      // -----------------------------------------------------------------
      // Hover tooltip + control toolbar
      // -----------------------------------------------------------------
      _setupTooltip(el, viewer, payload);
      _setupControls(el, viewer);

      return viewer;
    }

    // =========================================================================
    // Hover tooltip
    // =========================================================================

    function _setupTooltip(el, viewer, payload) {
      // Build a name → feature-spec lookup for tooltip content
      var lookup = {};
      var cgvData = payload && (payload.cgview || payload);
      var rawFeats = cgvData && cgvData.features;
      if (Array.isArray(rawFeats)) {
        rawFeats.forEach(function(f) {
          if (f.name) lookup[f.name] = f;
        });
      }

      // Tooltip div
      var tt = document.createElement('div');
      tt.style.cssText = [
        'position:absolute',
        'background:rgba(18,18,38,0.92)',
        'color:#f0f0f0',
        'padding:5px 10px',
        'border-radius:4px',
        'font:14px/1.6 monospace',
        'pointer-events:none',
        'opacity:0',
        'transition:opacity 0.1s',
        'z-index:200',
        'max-width:280px',
        'white-space:pre'
      ].join(';');
      el.style.position = 'relative';
      el.appendChild(tt);

      var canvas = el.querySelector('canvas');

      function _tooltipPos(e) {
        var er = el.getBoundingClientRect();
        var cr = canvas ? canvas.getBoundingClientRect() : er;
        var cx = (e.canvasX !== undefined ? cr.left + e.canvasX : e.clientX)
                  - er.left + 14;
        var cy = (e.canvasY !== undefined ? cr.top  + e.canvasY : e.clientY)
                  - er.top  - 10;
        return { x: cx, y: cy };
      }

      try {
        viewer.on('mousemove', function(e, data) {
          // CGView 1.8 passes hovered features in data.features or similar
          var hovered = null;
          if (data) {
            hovered = data.features
                   || (data.feature ? [data.feature] : null)
                   || (data.nearestFeature ? [data.nearestFeature] : null);
          }
          if (!hovered || !hovered.length) {
            tt.style.opacity = '0';
            return;
          }

          var f    = hovered[0];
          var spec = lookup[f.name] || {};
          var type = f.type   || spec.type   || '';
          var s    = f.start  || spec.start  || '';
          var e2   = f.stop   || spec.stop   || f.end || '';
          var str  = (f.strand === -1 || f.strand === '-1' || f.strand === '-')
                     ? '\u2212' : '+';

          tt.innerHTML =
            '<strong>' + _esc(f.name || '?') + '</strong>' +
            (type ? '  <span style="color:#aaa">[' + _esc(type) + ']</span>' : '') +
            '\n' + s + '..' + e2 + '  (' + str + ')';

          var p = _tooltipPos(e);
          tt.style.left    = p.x + 'px';
          tt.style.top     = p.y + 'px';
          tt.style.opacity = '1';
        });

        viewer.on('click', function() { tt.style.opacity = '0'; });
      } catch(e) {
        console.warn('[beta_cgview] could not attach hover events:', e);
      }

      el.addEventListener('mouseleave', function() { tt.style.opacity = '0'; });
    }

    // =========================================================================
    // Controls toolbar
    // =========================================================================

    function _setupControls(el, viewer) {
      var bar = document.createElement('div');
      bar.style.cssText = [
        'position:absolute',
        'top:10px',
        'left:10px',
        'display:flex',
        'flex-direction:column',
        'gap:3px',
        'z-index:100'
      ].join(';');

      var BTN = [
        'display:block',
        'width:34px',
        'height:34px',
        'border:none',
        'border-radius:5px',
        'background:rgba(30,30,50,0.78)',
        'color:#e8e8f0',
        'font:bold 16px/34px sans-serif',
        'text-align:center',
        'cursor:pointer',
        'padding:0',
        'user-select:none',
        'transition:background 0.12s'
      ].join(';');

      function makeBtn(label, tip, fn) {
        var b = document.createElement('button');
        b.textContent = label;
        b.setAttribute('title', tip);
        b.style.cssText = BTN;
        b.addEventListener('mouseenter', function() {
          b.style.background = 'rgba(30,30,50,0.97)';
        });
        b.addEventListener('mouseleave', function() {
          b.style.background = 'rgba(30,30,50,0.78)';
        });
        b.addEventListener('click', fn);
        return b;
      }

      function spacer() {
        var d = document.createElement('div');
        d.style.height = '4px';
        return d;
      }

      // ---- Zoom ----
      bar.appendChild(makeBtn('+', 'Zoom in', function() {
        try { viewer.zoomIn(); } catch(e) {}
      }));
      bar.appendChild(makeBtn('\u2212', 'Zoom out', function() {
        try { viewer.zoomOut(); } catch(e) {}
      }));
      bar.appendChild(spacer());

      // ---- Rotate ----
      bar.appendChild(makeBtn('\u25C4', 'Rotate left', function() {
        try { viewer.moveLeft(); } catch(e) {}
      }));
      bar.appendChild(makeBtn('\u25BA', 'Rotate right', function() {
        try { viewer.moveRight(); } catch(e) {}
      }));
      bar.appendChild(spacer());

      // ---- Reset ----
      bar.appendChild(makeBtn('\u21BA', 'Reset view', function() {
        try { viewer.reset(); } catch(e) {}
      }));

      // ---- Format toggle (circular ↔ linear) ----
      // Label shows what you'll switch TO: ○ = circular, ▬ = linear
      var isCircular = (viewer.format !== 'linear');
      var fmtBtn = makeBtn(
        isCircular ? '\u25AC' : '\u25EF',
        isCircular ? 'Switch to linear' : 'Switch to circular',
        null
      );
      fmtBtn.addEventListener('click', function() {
        try {
          var circ = (viewer.format !== 'linear');
          viewer.settings.update({ format: circ ? 'linear' : 'circular' });
          viewer.draw();
          fmtBtn.textContent = circ ? '\u25EF' : '\u25AC';
          fmtBtn.setAttribute('title', circ ? 'Switch to circular' : 'Switch to linear');
        } catch(e) {}
      });
      bar.appendChild(fmtBtn);
      bar.appendChild(spacer());

      // ---- Download PNG ----
      bar.appendChild(makeBtn('\u2B07', 'Download PNG (3000 px)', function() {
        try {
          var outH = 3000;
          var outW = Math.round(viewer.width / viewer.height * outH);
          viewer.io.downloadImage(outW, outH, 'cgview_map.png');
        } catch(e) {
          console.warn('[beta_cgview] download failed:', e);
        }
      }));

      el.appendChild(bar);
    }

    // =========================================================================
    // Utilities
    // =========================================================================

    function _esc(s) {
      return String(s)
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;');
    }

  }

});
