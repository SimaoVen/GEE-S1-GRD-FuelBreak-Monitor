/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var firebreakCollection = ee.FeatureCollection("projects/my-project-413614/assets/RPFGCs/ROI"),
    elevation = ee.Image("USGS/SRTMGL1_003");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
function powerToDb(img){
  return ee.Image(10).multiply(img.log10());
}

function dbToPower(img){
  return ee.Image(10).pow(img.divide(10));
}
 
// The RL speckle filter
function refinedLee(image) {
  var bandNames = image.bandNames();
  image = dbToPower(image);
  var result = ee.ImageCollection(bandNames.map(function(b){
    var img = image.select([b]);
    // img must be in natural units, i.e. not in dB!
    // Set up 3x3 kernels 
    var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
    var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);
    var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
    var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);
    // Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
    var sample_weights = ee.List([
      [0,0,0,0,0,0,0], 
      [0,1,0,1,0,1,0],
      [0,0,0,0,0,0,0], 
      [0,1,0,1,0,1,0], 
      [0,0,0,0,0,0,0], 
      [0,1,0,1,0,1,0],
      [0,0,0,0,0,0,0]]);
    var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);
    // Calculate mean and variance for the sampled windows and store as 9 bands
    var sample_mean = mean3.neighborhoodToBands(sample_kernel); 
    var sample_var = variance3.neighborhoodToBands(sample_kernel);
    // Determine the 4 gradients for the sampled windows
    var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
    gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
    gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
    gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());
    // And find the maximum gradient amongst gradient bands
    var max_gradient = gradients.reduce(ee.Reducer.max());
    // Create a mask for band pixels that are the maximum gradient
    var gradmask = gradients.eq(max_gradient);
    // duplicate gradmask bands: each gradient represents 2 directions
    gradmask = gradmask.addBands(gradmask);
    // Determine the 8 directions
    var directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
    directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
    directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
    directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
    // The next 4 are the not() of the previous 4
    directions = directions.addBands(directions.select(0).not().multiply(5));
    directions = directions.addBands(directions.select(1).not().multiply(6));
    directions = directions.addBands(directions.select(2).not().multiply(7));
    directions = directions.addBands(directions.select(3).not().multiply(8));
    // Mask all values that are not 1-8
    directions = directions.updateMask(gradmask);
    // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
    directions = directions.reduce(ee.Reducer.sum());  
    var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));
    // Calculate localNoiseVariance
    var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);
    // Set up the 7*7 kernels for directional statistics
    var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));
    var diag_weights = ee.List([
      [1,0,0,0,0,0,0], 
      [1,1,0,0,0,0,0], 
      [1,1,1,0,0,0,0], 
      [1,1,1,1,0,0,0], 
      [1,1,1,1,1,0,0], 
      [1,1,1,1,1,1,0], 
      [1,1,1,1,1,1,1]]);
    var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
    var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);
    // Create stacks for mean and variance using the original kernels. Mask with relevant direction.
    var dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
    var dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));
    // and add the bands for rotated kernels
    for (var i=1; i<4; i++) {
      dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
      dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
      dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
      dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    }
    // "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
    dir_mean = dir_mean.reduce(ee.Reducer.sum());
    dir_var = dir_var.reduce(ee.Reducer.sum());
    // A finally generate the filtered value
    var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));
    b = varX.divide(dir_var);
    return dir_mean.add(b.multiply(img.subtract(dir_mean)))
      .arrayProject([0])
      // Get a multi-band image bands.
      .arrayFlatten([['sum']])
      .float();
  })).toBands().rename(bandNames);
  return powerToDb(result);
}

var ninetyRad = ee.Image.constant(90).multiply(Math.PI/180);
// Volumetric Model Hoekman 1990
function volumetricModel(theta_iRad, alpha_rRad){
  var nominator = (ninetyRad.subtract(theta_iRad).add(alpha_rRad)).tan();
  var denominator = (ninetyRad.subtract(theta_iRad)).tan();
  return nominator.divide(denominator);
}
  
// Calculate masks
function _masking(alpha_rRad, theta_iRad, proj){
  // layover, where slope > radar viewing angle
  var layover = alpha_rRad.lt(theta_iRad).rename('layover');
  // shadow
  var shadow = alpha_rRad.gt(ee.Image.constant(-1).multiply(ninetyRad.subtract(theta_iRad))).rename('shadow');
  // combine layover and shadow
  var mask = layover.and(shadow);
  return mask.rename('no_data_mask');
}

function slope_correct(image){
  // get image geometry and projection
  var geom = image.geometry();
  var proj = image.select(1).projection();
  // get look direction angle
  var heading = (ee.Terrain.aspect(image.select('angle')).reduceRegion(ee.Reducer.mean(), geom, 1000).get('aspect'));
  // Sigma0 to Power of input image
  var sigma0Pow = dbToPower(image);
  
  // Radar geometry
  var theta_iRad = image.select('angle').multiply(Math.PI/180).clip(geom);
  var phi_iRad = ee.Image.constant(heading).multiply(Math.PI/180);
  
  // Terrain geometry
  var alpha_sRad = ee.Terrain.slope(elevation).select('slope').multiply(Math.PI/180).setDefaultProjection(proj).clip(geom);
  var phi_sRad = ee.Terrain.aspect(elevation).select('aspect').multiply(Math.PI/180).setDefaultProjection(proj).clip(geom);

  // Model geometry
  // reduce to 3 angle
  var phi_rRad = phi_iRad.subtract(phi_sRad);
  // slope steepness in range
  var alpha_rRad = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan();
  // slope steepness in azimuth
  var alpha_azRad = (alpha_sRad.tan().multiply(phi_rRad.sin())).atan();
  // Gamma_nought
  var gamma0 = sigma0Pow .divide(theta_iRad.cos());

  // Apply the model
  var corrModel = volumetricModel(theta_iRad, alpha_rRad);

  // apply model to derive gamma0_flat
  var gamma0_flat = gamma0.divide(corrModel);
  // transform to dB-scale
  var gamma0_flatDB = powerToDb(gamma0_flat).select(['VV', 'VH']);
  // get Layover/Shadow mask
  var mask = _masking(alpha_rRad, theta_iRad, proj);
  // return gamma_flat plus mask
  return ee.Image(gamma0_flatDB.addBands(mask).copyProperties(image));
}

function getCollection(year, geometry){
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate = startDate.advance(1, 'year');
  var collection = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .filter(ee.Filter.date(startDate, endDate))
    .filterBounds(geometry);
  return collection;
}

function correctImages(year, collection, zoneID) {
  var bands = collection.first().bandNames();
  var sequenceEnd = 11;
  var skipMode = 'month';
  var increment = 1;
  return ee.ImageCollection(ee.List.sequence(0, sequenceEnd).map(function(time) {
    var advanceAmount = ee.Number(time).multiply(increment);
    var startDate = ee.Date.fromYMD(year, 1, 1).advance(advanceAmount, skipMode);
    var endDate = startDate.advance(increment, skipMode);
    var corrected = collection.filterDate(startDate, endDate).map(function(image) {
      var slope_corrected = slope_correct(image);
      var lee_filtered = refinedLee(slope_corrected);
      return lee_filtered;
    })
      .reduce(ee.Reducer.median())
      .rename(bands);
    return corrected.set('cod_troco', zoneID).set('system:time_start', startDate.millis());
  }));
}

function createRatioBand(collection) {
  return collection.map(function(img) {
    var ratio = img.select('VV').divide(img.select('VH')).rename('ratio');
    return img.addBands(ratio);
  });
}

function makeCollection(year, geometry, zoneID) {
  var raw_collection = getCollection(year, geometry);
  var corrected_collection = createRatioBand(correctImages(year, raw_collection, zoneID));
  var colList = corrected_collection.toList(corrected_collection.size());
  var n = colList.size().getInfo();
  for (var i = 0; i < n; i++) {
    var img = ee.Image(colList.get(i));
    var id = 'Images/Monthly/' + zoneID + '_' + ee.Date(img.get('system:time_start')).format('yyyyMM').getInfo();
    print(id);
    Export.image.toAsset({
      image:img,
      description: i.toString(),
      assetId: id,
      region: geometry,
      scale: 20});
  }
}

Map.add(ui.Map.Layer(firebreakCollection, {}, 'Firebreak Collection'));
var panel = ui.Panel({ 
  style: { 
    position: 'top-center', 
    padding: '10px' 
  }
}).setLayout(ui.Panel.Layout.flow('horizontal'));

var selectedYear;
var yearSelector = ui.Select({ 
  items: ['2020', '2021', '2022', '2023', '2024'], 
  placeholder: 'Select a Year', 
  onChange: function(year) {selectedYear = ee.Number.parse(year);} 
});

var selectedFB;
var uniqueValuesFB = firebreakCollection.aggregate_array('cod_troco').distinct().sort();
uniqueValuesFB.evaluate(function(values) {
  var fbSelector = ui.Select({ 
    items: values, 
    placeholder: 'Select a FB', 
    onChange: function(value) {
      selectedFB = value;
    }
  });
  panel.add(fbSelector);
});

var button = ui.Button({ 
  label: 'Process Images', 
  onClick: function() {
    if(!selectedYear || !selectedFB)
      return;
    makeCollection(selectedYear, firebreakCollection.filter(ee.Filter.eq('cod_troco', selectedFB)).geometry().bounds(), selectedFB);
  }
});

panel.add(button);
panel.add(yearSelector);
Map.add(panel);