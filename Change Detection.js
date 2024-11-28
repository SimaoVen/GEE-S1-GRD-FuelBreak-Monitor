/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var sarVisParams = {"min":[-15,-20,1],"max":[0,-5,15],"bands":["VV","VH","ratio"]},
    fuelbreakCollection = ee.FeatureCollection("projects/my-project-413614/assets/RPFGCs/FB");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
function get_collection(year, zone) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate = startDate.advance(1, 'year');
  var imageList = ee.data.listAssets('projects/my-project-413614/assets/Images/Monthly').assets.map(function(d) { return d.name });
  return ee.ImageCollection(imageList)
    .filter(ee.Filter.eq('cod_troco', zone))
    .filter(ee.Filter.date(startDate, endDate))
    .sort('system:time_start', true);
}

var binary_classifier = ee.Classifier.load('projects/my-project-413614/assets/binary_classifier');
function get_binary(collection) {
  var classified_collection = collection.map(function(image) {
    return image.classify(binary_classifier)
      .focalMin({kernel: ee.Kernel.circle({radius: 10, units: 'meters'}), iterations: 5})
      .focalMax({kernel: ee.Kernel.circle({radius: 10, units: 'meters'}), iterations: 5})
      .copyProperties(image, ['system:time_start']);
  });
  return classified_collection;
}

function get_change(collection) {
  var collectionList = collection.toList(collection.size());
  var indices =  ee.List.sequence(1, collection.size().subtract(1));
  return ee.ImageCollection(indices.map(function(index) {
    var image1 = ee.Image(collectionList.get(ee.Number(index).subtract(1)));
    var image2 = ee.Image(collectionList.get(index));
    return image1.subtract(image2).copyProperties(image2, ['system:time_start']);
  }));
}

function get_change_evolution(collection) {
  var collectionSize = collection.size();
  var collectionList = collection.toList(collectionSize);
  var firstOccurrence = ee.Image.constant(0);
  var sequence = ee.List.sequence(0, collectionSize.subtract(1));
  var firstOccurrenceImage = ee.Image(sequence.iterate(function(index, image) {
    index = ee.Number(index);
    image = ee.Image(image);
    // Cast the current image in the collection to check for the value 1
    var currentImage = ee.Image(collectionList.get(index));
    // Create a condition mask where the current image has value 1 and the firstOccurrence image is still 0
    var condition = currentImage.eq(1).and(image.eq(0));
    // Update the firstOccurrence image with the index value where the condition is true
    return image.where(condition, index.add(1));
  }, firstOccurrence));
  // Mask the firstOccurrenceImage to only show the index of the first occurrence of value 1
  var maskedImage = firstOccurrenceImage.updateMask(firstOccurrenceImage.gt(0));
  return maskedImage;
}


// Map
var fuelBreakCollectionLayer = ui.Map.Layer(fuelbreakCollection, {}, 'Fuel Break Collection');
var fuelBreakLayer;
var compositeLayer1;
var compositeLayer2;
var changeLayer;
var changeEvoLayer;
var mapPanel = ui.Map();
var mapLayers = mapPanel.layers();
mapLayers.add(fuelBreakCollectionLayer);

// Panel Setup
var userPanel = ui.Panel({style: {width: '35%'}});

var intro = ui.Panel([
  ui.Label({
    value: 'SAR Portugal Fuel Break Change Detection | Time Series Inspector',
    style: {fontSize: '20px', fontWeight: 'bold'}
  }),
  ui.Label('- Select a fuel break from the map, the year and region to generate the time series charts.')
]);
userPanel.add(intro);

var fuelBreakInfo = ui.Label('No Fuel Break Selected');
userPanel.add(ui.Panel([ui.Label('Fuel Break Info: '), fuelBreakInfo], ui.Panel.Layout.flow('horizontal')));
var selectedFuelBreak;
function updateSelectedFuelBreak(coords) {
  // Create a point geometry from the clicked coordinates
  var point = ee.Geometry.Point(coords.lon, coords.lat);
  // Filter the Fuel Break FeatureCollection
  selectedFuelBreak = fuelbreakCollection.filterBounds(point);
  // Remove the last added layer
  if(fuelBreakLayer)
    mapLayers.remove(fuelBreakLayer);
  // Add the selected feature to the map with the defined style and update roi
  fuelBreakLayer = ui.Map.Layer(selectedFuelBreak, {color: 'FF5F00'}, 'Selected Fuel Break');
  mapLayers.add(fuelBreakLayer);
  var id = selectedFuelBreak.first().get('id');
  var cod_troco = selectedFuelBreak.first().get('cod_troco');
  id.evaluate(function(result) {
    cod_troco.evaluate(function(result1) {
      if(result || result1) {
        fuelBreakInfo.setValue('id: ' + result + ' | cod_troco: ' + result1);
      }
      else {
        fuelBreakInfo.setValue('No Fuel Break Selected');
      }
    });
  });
}

var collection; // S1-GRD Processed Images
var bin_collection; // Binary Classified Images
var change_collection;  // Change Maps
var change_evo; // Yearly Change Map
function updateCollections(year, zone) {
  if(!year || !zone)
    return;
  collection = get_collection(year, zone);
  bin_collection= get_binary(collection);
  change_collection = get_change(bin_collection);
  change_evo = get_change_evolution(change_collection);
}

var years = ['2021', '2022', '2023'];
var select_year = ui.Select({
  items: years,
  placeholder: 'Select Year',
  onChange: function(new_year) {
    updateCollections(ee.Number.parse(new_year), select_zone.getValue());
  }
});

var zones = ['BOTICAS', 'SEIA'];
var select_zone = ui.Select({
  items: zones,
  placeholder: 'Select ROI',
  onChange: function(zone) {
    updateCollections(ee.Number.parse(select_year.getValue()), zone);
  }
});

var chartPanel = ui.Panel();

function get_date_string(date) {
  return date.format('MMMM yyyy');
}

function sar_chart() {
  var chart_title = 'Mean Backscatter Time Series';
  var chart = ui.Chart.image.series({
    imageCollection: collection.select(['VV', 'VH', 'ratio']),
    region: selectedFuelBreak.geometry(),
    reducer: ee.Reducer.mean(),
    scale: 20,
    xProperty: 'system:time_start'
  });
  // Rename the series bands
  chart.setSeriesNames(['VV', 'VH', 'VV/VH']);
  chart.setOptions({
    title: chart_title,
    hAxis: {
      title: 'Date',
      titleTextStyle: {italic: false, bold: true},
      format: 'MM-yyyy'
    },
    vAxis: {
      title: 'Backscatter (dB)',
      titleTextStyle: {italic: false, bold: true}
    },
    lineWidth: 5,
    colors: ['#ff8000', '#8000ff', '#ad5050'],
    curveType: 'function'
  });
  return chart;
}

function change_chart() {
  var classNames = ee.Dictionary.fromLists(['-1', '0', '1'], ['Change to Vegetation', 'No Change', 'Change to Bare Ground']);
  var classColors =  ee.Dictionary.fromLists(['-1', '0', '1'], ['#98ff00', '#000000', '#c57704']);
  var chart_title = 'Change Classification';
  var classInfoList = change_collection.map(function(image){
    var areas = ee.Image.pixelArea().divide(1e4).addBands(image).reduceRegion({
      reducer: ee.Reducer.sum().group({
        groupField: 1,
        groupName: 'classification'
      }),
      geometry: selectedFuelBreak.geometry(),
      scale: 20
    });
    var classInfo = ee.List(areas.get('groups')).map(function(item) {
      var areaDict = ee.Dictionary(item);
      var classNumber = areaDict.getNumber('classification').format();
      var className = classNames.get(classNumber);
      var classArea = areaDict.getNumber('sum');
      var classColor = classColors.get(classNumber);
      var time = image.get('system:time_start');
      return ee.Feature(selectedFuelBreak.geometry(), {
        'class': classNumber,
        'className': className,
        'area': classArea,
        'color': classColor,
        'system:time_start': time
      });
    });
    return ee.FeatureCollection(classInfo);
  }).flatten();
  var classInfoFeatureCollection = ee.FeatureCollection(classInfoList.toList(classInfoList.size()).map(function(element) {
    return ee.Feature(element);
  }));
  var chart = ui.Chart.feature.groups({
    features: classInfoFeatureCollection,
    xProperty: 'system:time_start',
    yProperty: 'area',
    seriesProperty: 'className'
  }).setChartType('LineChart');
  classInfoFeatureCollection.aggregate_array('color').evaluate(function(colorList){
    chart.setOptions({
      title: chart_title + ' Time Series',
      hAxis: {
        title: 'Date',
        titleTextStyle: {italic: false, bold: true},
        format: 'MM-yyyy'
      },
      vAxis: {
        title: 'Area (ha)',
        titleTextStyle: {italic: false, bold: true}
      },
      colors: colorList
    });
  });
  chart.onClick(function(xValue){
    if (!xValue) return;  // Selection was cleared.
    var col = collection.sort('system:time_start', false);
    var SAR_image = ee.Image(col.filter(ee.Filter.eq('system:time_start', xValue)).first());
    var SAR_previoius = ee.Image(col.filter(ee.Filter.lt('system:time_start', xValue)).first());
    var Change_image = ee.Image(change_collection.filter(ee.Filter.eq('system:time_start', xValue)).first());
    var date1 = get_date_string(ee.Date(SAR_image.get('system:time_start')));
    var date2 = get_date_string(ee.Date(SAR_previoius.get('system:time_start')));
    if(changeLayer)
      mapLayers.remove(changeLayer);
    changeLayer = ui.Map.Layer(Change_image,{min:-1, max:1, palette: ['#98ff00', '#000000', '#c57704']}, chart_title + ' for ' + date1.getInfo(), false);
    mapLayers.add(changeLayer);
    if(compositeLayer1)
      mapLayers.remove(compositeLayer1);
    compositeLayer1 = ui.Map.Layer(SAR_image, sarVisParams, 'SAR Composite for ' + date1.getInfo(), false);
    mapLayers.add(compositeLayer1);
    if(compositeLayer2)
      mapLayers.remove(compositeLayer2);
    compositeLayer2 = ui.Map.Layer(SAR_previoius, sarVisParams, 'SAR Composite for ' + date2.getInfo(), false);
    mapLayers.add(compositeLayer2);
    //if(changeEvoLayer)
    //  mapLayers.remove(changeEvoLayer);
    changeEvoLayer = ui.Map.Layer(change_evo.clip(selectedFuelBreak.geometry()),{min: 1,max: 12,palette: ['#484cbc', '#4884fc', '#30bcec', '#20e4b4', '#60fc74', '#a8fc3c', '#e0e434', '#ffbc3c', '#ff7c24', '#e8440c', '#c01c04', '#800404']}, 'Change Evo', false);
    mapLayers.add(changeEvoLayer);
    Export.image.toDrive({
      image: change_evo,
      description: 'image_export',
      folder: 'GEE',
      region: fuelbreakCollection.filter(ee.Filter.eq('cod_troco', select_zone.getValue())).geometry().bounds(),
      scale: 20,
      crs: 'EPSG:4326'
    });
  });
  return chart;
}

function lcc_chart() {
  var classNames = ee.Dictionary.fromLists(['0', '1'], ['Bare Ground', 'Vegetation']);
  var classColors =  ee.Dictionary.fromLists(['0', '1'], ['#c57704', '#98ff00']);
  var chart_title = 'Land Cover Area';
  var classInfoList = bin_collection.map(function(image){
    var areas = ee.Image.pixelArea().divide(1e4).addBands(image).reduceRegion({
      reducer: ee.Reducer.sum().group({
        groupField: 1,
        groupName: 'classification'
      }),
      geometry: selectedFuelBreak.geometry(),
      scale: 20
    });
    var classInfo = ee.List(areas.get('groups')).map(function(item) {
      var areaDict = ee.Dictionary(item);
      var classNumber = areaDict.getNumber('classification').format();
      var className = classNames.get(classNumber);
      var classArea = areaDict.getNumber('sum');
      var classColor = classColors.get(classNumber);
      var time = image.get('system:time_start');
      return ee.Feature(selectedFuelBreak.geometry(), {
        'class': classNumber,
        'className': className,
        'area': classArea,
        'color': classColor,
        'system:time_start': time
      });
    });
    return ee.FeatureCollection(classInfo);
  }).flatten();
  var classInfoFeatureCollection = ee.FeatureCollection(classInfoList.toList(classInfoList.size()).map(function(element) {
    return ee.Feature(element);
  }));
  var chart = ui.Chart.feature.groups({
    features: classInfoFeatureCollection,
    xProperty: 'system:time_start',
    yProperty: 'area',
    seriesProperty: 'className'
  }).setChartType('LineChart');
  classInfoFeatureCollection.aggregate_array('color').evaluate(function(colorList){
    chart.setOptions({
      title: chart_title + ' Time Series',
      hAxis: {
        title: 'Date',
        titleTextStyle: {italic: false, bold: true},
        format: 'MM-yyyy'
      },
      vAxis: {
        title: 'Area (ha)',
        titleTextStyle: {italic: false, bold: true}
      },
      colors: colorList
    });
  });
  return chart;
}

// Create a button to generate a chart.
var create_chart = ui.Button({
  label: 'Generate Chart',
  onClick: function() {
    chartPanel.clear();
    if(fuelBreakInfo.getValue() == 'No Fuel Break Selected') {
      chartPanel.add(ui.Label('Please select a Fuel Break'));
      return;
    }
    chartPanel.add(sar_chart());
    chartPanel.add(lcc_chart());
    chartPanel.add(change_chart());
  }
});

// Create a button that clears the panel when clicked.
var reset = ui.Button({
  label: 'Reset Map Layers',
  onClick: function() {
    if(compositeLayer1)
      mapLayers.remove(compositeLayer1);
    if(compositeLayer1)
      mapLayers.remove(compositeLayer1);
    if(changeLayer)
      mapLayers.remove(changeLayer);
    //if(changeEvoLayer)
    //  mapLayers.remove(changeEvoLayer);
  }
});

userPanel.add(ui.Panel([select_year, select_zone, create_chart, reset], ui.Panel.Layout.flow('horizontal')));
userPanel.add(chartPanel);
// Center map to Portugal.
mapPanel.centerObject(ee.Geometry.Point(-8.1609, 39.1163), 7);
// When the map is clicked try to update the selected fuel break.
mapPanel.onClick(updateSelectedFuelBreak);
// Replace the root with a SplitPanel that contains the inspector and map.
ui.root.clear();
ui.root.add(ui.SplitPanel(userPanel, mapPanel));