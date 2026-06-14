var Module = {
  onRuntimeInitialized: function () {
    var btn = document.getElementById('price-btn');
    btn.textContent = 'Price It';
    btn.disabled = false;
    btn.addEventListener('click', onPriceClick);
  }
};

var spatialResSlider = document.getElementById('spatial-res');
var temporalResSlider = document.getElementById('temporal-res');
spatialResSlider.addEventListener('input', function () {
  document.getElementById('spatial-res-val').textContent = this.value;
});
temporalResSlider.addEventListener('input', function () {
  document.getElementById('temporal-res-val').textContent = this.value;
});

function validateParams(strike, expiration, spot, sigma, rate, dividendYield, spatialRes, timeSteps) {
  if (sigma <= 0) return 'Volatility must be positive.';
  if (rate < 0) return 'Risk-free rate cannot be negative.';
  if (dividendYield < 0) return 'Dividend yield cannot be negative.';
  if (expiration <= 0) return 'Expiration must be positive.';
  if (spot <= 0) return 'Underlying price must be positive.';
  if (strike <= 0) return 'Strike price must be positive.';
  if (spatialRes < 5) return 'Spatial nodes must be at least 5.';
  if (timeSteps < 2) return 'Time steps must be at least 2.';
  return null;
}

function onPriceClick() {
  var btn = document.getElementById('price-btn');
  var loading = document.getElementById('loading');
  var errorDiv = document.getElementById('error-display');
  var priceDisplay = document.getElementById('price-display');

  var optionType, strike, expiration, spot, sigma, rate, dividendYield, spatialRes, timeSteps, pricingModel;
  try {
    pricingModel = parseInt(document.getElementById('pricing-model').value);
    optionType = parseInt(document.getElementById('option-type').value);
    strike = parseFloat(document.getElementById('strike').value);
    expiration = parseFloat(document.getElementById('expiration').value) / 365;
    spot = parseFloat(document.getElementById('spot').value);
    sigma = parseFloat(document.getElementById('sigma').value) / 100;
    rate = parseFloat(document.getElementById('rate').value) / 100;
    dividendYield = parseFloat(document.getElementById('dividend-yield').value) / 100;
    spatialRes = parseInt(spatialResSlider.value);
    timeSteps = parseInt(temporalResSlider.value);
  } catch (e) {
    errorDiv.textContent = e.message || 'Failed to read parameters.';
    errorDiv.classList.remove('hidden');
    return;
  }

  var error = validateParams(strike, expiration, spot, sigma, rate, dividendYield, spatialRes, timeSteps);
  if (error) {
    errorDiv.textContent = error;
    errorDiv.classList.remove('hidden');
    return;
  }
  errorDiv.classList.add('hidden');

  btn.disabled = true;
  loading.classList.remove('hidden');
  priceDisplay.textContent = 'Option Price: computing...';

  setTimeout(function () {
    try {
      var handle = Module._pricer_create();
      if (!handle) {
        throw new Error('Failed to create pricer (out of memory).');
      }

      var surfSize = timeSteps * spatialRes;
      var surfPtr = Module._malloc(surfSize * 8);
      var sGridPtr = Module._malloc(spatialRes * 8);
      var tGridPtr = Module._malloc(timeSteps * 8);

      if (!surfPtr || !sGridPtr || !tGridPtr) {
        if (surfPtr) Module._free(surfPtr);
        if (sGridPtr) Module._free(sGridPtr);
        if (tGridPtr) Module._free(tGridPtr);
        Module._pricer_destroy(handle);
        throw new Error('Failed to allocate WASM memory for grids.');
      }

      Module._pricer_surface(handle, pricingModel, optionType, strike, expiration,
                             spot, dividendYield, sigma, rate, timeSteps, spatialRes,
                             surfPtr, sGridPtr, tGridPtr);

      var price = Module._pricer_price(handle, pricingModel, optionType, strike, expiration,
                                       spot, dividendYield, sigma, rate, timeSteps, spatialRes);

      var sGrid = Array.from(new Float64Array(Module.HEAPU8.buffer, sGridPtr, spatialRes));
      var tGrid = Array.from(new Float64Array(Module.HEAPU8.buffer, tGridPtr, timeSteps));

      var surfaceFlat = new Float64Array(Module.HEAPU8.buffer, surfPtr, surfSize);
      var zData = [];
      for (var i = 0; i < timeSteps; i++) {
        zData.push(Array.from(surfaceFlat.subarray(i * spatialRes, (i + 1) * spatialRes)));
      }

      Module._free(surfPtr);
      Module._free(sGridPtr);
      Module._free(tGridPtr);
      Module._pricer_destroy(handle);

      priceDisplay.textContent = 'Option Price: ' + price.toFixed(4);

      Plotly.newPlot('surface-plot', [{
        type: 'surface',
        x: sGrid,
        y: tGrid,
        z: zData,
        // colorscale: 'Viridis',
        colorscale: 'Cividis',
        contours: {
            z: { show: true, usecolormap: true, highlightcolor: 'rgba(255, 136, 0, 1)', project: { z: true } },
          x: { show: true, usecolormap: true, highlightcolor: 'rgba(179, 255, 0, 1)', project: { x: true } }
        }
      }], {
        title: (pricingModel === 0 ? 'European' : 'American') + ' Option Price Surface',
        scene: {
          xaxis: { title: 'Stock Price', color: '#e0e0e0', gridcolor: '#2a2a4a' },
          yaxis: { title: 'Time to Maturity (years)', color: '#e0e0e0', gridcolor: '#2a2a4a' },
          zaxis: { title: 'Option Price', color: '#e0e0e0', gridcolor: '#2a2a4a' }
        },
        margin: { l: 0, r: 0, b: 0, t: 50 },
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        font: { color: '#e0e0e0', family: '-apple-system, BlinkMacSystemFont, Segoe UI, Roboto, sans-serif' }
      }, { responsive: true });

    } catch (e) {
      errorDiv.textContent = e.message || 'An unexpected error occurred.';
      errorDiv.classList.remove('hidden');
      priceDisplay.textContent = 'Option Price: &mdash;';
    }

    loading.classList.add('hidden');
    btn.disabled = false;
  }, 50);
}
