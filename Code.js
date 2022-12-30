// Display the result.
// Map.centerObject(image, 10);
// Map.addLayer(maxValue, {max: 13000}, 'Maximum value image');


// NDVI
var rgb_vis = {min: 0, max: 0.3, bands: ['B4', 'B3', 'B2']};
var filtered = imageCollection2.filterDate('2018-06-01', '2018-07-30')
  .filterBounds(roi);
var image = ee.Image(filtered.first());
var ndvi = image.normalizedDifference(['B5', 'B4']);
print(ui.Chart.image.series(imageCollection2, AOI_subset,ee.Reducer.mean(), 30));
Map.addLayer(image, rgb_vis, 'RGB');
Map.addLayer(ndvi, {min: 0, max: 1}, 'NDVI');
var vizParams = {
  bands: ['B5', 'B4', 'B3'],
  min: 926,
  max: 6000,
  gamma: [0.95, 1.1, 1]
};
// Map.addLayer(image,vizParams,"Colors");

// Load an image and select some bands of interest.
// var image = ee.Image('LANDSAT/LC08/C01/T1/LC08_044034_20140318')
//     .select(['B4', 'B3', 'B2']);
    
var dataset = ee.FeatureCollection('USGS/WBD/2017/HUC12');
  // Map.addLayer(dataset);
var fc= dataset.filter(ee.Filter.eq("huc12","101900071002"));
Map.addLayer(fc);
var Landcover = ee.ImageCollection('USGS/NLCD')
                      .map(function(image){
                      return image.clipToCollection(fc).select("landcover");
                      });
                      
// var NLCD_Fossil = Landcover.select('landcover');

// Load an image collection, filtered so it's not too much data.
var collection = ee.ImageCollection('LANDSAT/LT05/C01/T1')
  .filterDate('2008-01-01', '2008-12-31')
  .filter(ee.Filter.eq('WRS_PATH', 44))
  .filter(ee.Filter.eq('WRS_ROW', 34));

// Compute the median in each band, each pixel.
// Band names are B1_median, B2_median, etc.
var medianlc = Landcover.reduce(ee.Reducer.median());

// The output is an Image.  Add it to the map.
var vis_param = {bands: ['landcover_median'], gamma: 1.6};
Map.setCenter(-122.3355, 37.7924, 9);
Map.addLayer(medianlc, vis_param,"Landcover NDVI median");
    

var ndvinlcd = ndvi.addBands(medianlc);
// Reduce the image to get a one-band maximum value image.
 var means = ndvinlcd.reduceRegion({
  reducer: ee.Reducer.mean().group({
    groupField: 1,
    groupName: 'code',
  }),
  // geometry: region.geometry(),
  scale: 1000,
  maxPixels: 1e8
});

// Print the resultant Dictionary.
print(means);
 
 
 
Map.setCenter(-105.008849,40.490646,12);
// Map.addLayer(NLCD_Fossil,{},"Fossil NLCD");
Map.addLayer(Landcover.first(),{},"NDVINLCD");
// Map.addLayer(Landcover.first(),{},"Landcover");
// Map.addLayer(medianlc,{},"NdVI");

//////////////////////////////////////////////////////////////
// Asset List
//////////////////////////////////////////////////////////////
// Outline for AOI (subset); loaded as asset on user account;
// A general polygon will work for testing purposes


// Image collection for Landsat 8 TOA
var l8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA');
// Filter by AOI
var spatialFiltered = l8.filterBounds(AOI_subset);
// Filter by date
var temporalFiltered = spatialFiltered.filterDate('2017-07-31', '2018-07-31');
// Filter by cloud cover
var sorted = temporalFiltered.filter(ee.Filter.lt('CLOUD_COVER', 10));

// Get the median over time, in each band, in each pixel.
var median = temporalFiltered.median();

// Create a function to calculate NDWI per image
function calcNDWI(image){
  return image.normalizedDifference(['B3', 'B5']).rename('NDWI');
}

//////////////////////////////////////////////////////////////
// Otsu Segmentation
//////////////////////////////////////////////////////////////

// Use Otsu (1979) method for image segmentation to create water-land mask
// Based on code from Nicholas Clinton (Google): https://code.earthengine.google.com/6f3abb73b6642198485e36000024827a
// Return the DN that maximizes interclass variance in B5 (in the region).
var otsu = function(histogram) {
  var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
  var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);
  
  var indices = ee.List.sequence(1, size);
  
  // Compute between sum of squares, where each mean partitions the data.
  var bss = indices.map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMeans = means.slice(0, 0, i);
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
           bCount.multiply(bMean.subtract(mean).pow(2)));
  });
  
  // Return the mean value corresponding to the maximum BSS.
  return means.sort(bss).get([-1]);
};

//////////////////////////////////////////////////////////////
// Map over the entire collection
//////////////////////////////////////////////////////////////
var mappedCollection = sorted.map(function(image){
  var NDWI = calcNDWI(image);
  var histogram = ee.Dictionary(NDWI.reduceRegion({
    reducer: ee.Reducer.histogram(100),
    geometry: AOI_subset, 
    scale: 30,
    bestEffort: true
  }).get('NDWI'));
  var threshold = otsu(histogram);
  var waterMask = NDWI.lt(ee.Number(threshold)).rename('waterMask');
  return image.addBands([NDWI, waterMask]).set('OtsuThresh', threshold);
});


//////////////////////////////////////////////////////////////
// Do Otsu threshold for the median image
//////////////////////////////////////////////////////////////

var NDWImedian = calcNDWI(median);
var histogram = ee.Dictionary(NDWImedian.reduceRegion({
  reducer: ee.Reducer.histogram(100),
  geometry: AOI_subset, 
  scale: 30,
  bestEffort: true
}).get('NDWI'));
var threshold = otsu(histogram);
var waterMask = NDWImedian.lt(ee.Number(threshold)).rename('waterMask');
var median = median.addBands([NDWImedian, waterMask]).set('OtsuThresh', threshold);

//////////////////////////////////////////////////////////////
// Visualization and Printing
//////////////////////////////////////////////////////////////
// Print
print('size full collection', sorted.size());
print('The full collection', mappedCollection);
print('Otsu thresholds', mappedCollection.aggregate_array('OtsuThresh'));
print("Median NDWI of Fossil Creek Reservoir",ui.Chart.image.histogram(median.select("NDWI"),AOI_subset,30));

// Make a variable of visualization parameters.
var visParams = {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.35};
var ndwiViz = {min: 0.01, max: 0.75, palette: ['2ca25f', '0000FF'], bands: 'NDWI'};
/////////////////////
////////////////////
//////////////////
//////////////////////
////////////////////

// Load NLCD landcover map and clip to watershed
var nlcd = ee.ImageCollection('USGS/NLCD');

nlcd = nlcd.toList(nlcd.size());
print('dataset: ', nlcd);
nlcd = ee.ImageCollection([nlcd.get(1), nlcd.get(5),nlcd.get(6),nlcd.get(7),
        nlcd.get(8),nlcd.get(12),nlcd.get(13)]);
// Select the 'landcover' band
var wland = nlcd.select('landcover');
print('Landcover: ',wland);


Map.addLayer(median.updateMask(median.select('waterMask')), visParams, 'True color median water masked');
Map.addLayer(mappedCollection.first().updateMask(mappedCollection.first().select('waterMask')), visParams, 'True color first image water masked');
Map.addLayer(median, ndwiViz, 'NDWI median', 0);
Map.addLayer(mappedCollection.first(), ndwiViz, 'NDWI first image', 0);

