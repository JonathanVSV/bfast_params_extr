// -------Imports section
var l5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR"),
    l7 = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR"),
    l8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR");

// --------Constants definition--------------
var folder1 = 'Bfast_Ayuquila',
    // Period of images that are going to be used
    dateIniIm = '1990-08-01',
    dateEndIm = '2019-01-01',
    // Path row
    path = 29,
    row = 46,
    // Maximum cloud cover for images in image collection
    maxCCL = 100,
    // Resolution in m for export
    resolution = 30,
    // Region of interest and rectangle covering all the roi
    roi = disturbAreas,
    rectangle = geometry;

//--------- Functions-------------
//Rename function  
  var renameFunc = function(image) {
    return image
    .rename('B','G','R','NIR','SWIR1','SWIR2');
  };

//Filter images function
var filterImCol = function (imageCol,sensor,dateInicial,dateFinal){
 
  var bandas = ee.Algorithms.If(sensor == 8,
    ['B2', 'B3', 'B4', 'B5', 'B6', 'B7'],
    ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']);
    
  // Clear observations for Land, water is ignored
  var CLEAR = ee.Algorithms.If(sensor == 8, 322, 66);

  var maskCloudsWrap = function(CLEAR){
    var maskClouds = function(image) {
      var fMask = image.select('pixel_qa');
        
      var cloudMask = fMask.eq(ee.Number(CLEAR));
      
      return image.updateMask(cloudMask);
    };
    return maskClouds;
  };
  
  // Apply filters and mask clouds
  var temp = imageCol.filterDate(dateInicial,dateFinal)
   //Area of interest filter
   .filterBounds(rectangle)
   //Cloud cover of the image filter
   .filter(ee.Filter.lte('CLOUD_COVER_LAND',maxCCL))
   .filter(ee.Filter.eq('WRS_PATH', path))
   .filter(ee.Filter.eq('WRS_ROW', row))
   //sd´í msdkClouds ya hace cosas diferentes dependiendeo del sensor
   .map(maskCloudsWrap(CLEAR))
   .select(bandas);
 
    //Rename bands from B1, B2, etc, to R, NIR, SWIR, etc
    temp = temp.map(renameFunc);
    temp = temp.map(function(image){
                      return image.normalizedDifference(['NIR','R']).rename('ndvi')});

    return temp;
  };

// ---------------Use functions-----------------------------
// Filter landsat images and mask clouds
var l5Filtered = filterImCol(l5, 5, dateIniIm, dateEndIm);
var l7Filtered = filterImCol(l7, 7, dateIniIm, dateEndIm);
var l8Filtered = filterImCol(l8, 8, dateIniIm, dateEndIm);

// Merge 3 collection
var lAll = l5Filtered.merge(l7Filtered).merge(l8Filtered);

// Transform images into int16 to reduce size.
lAll = lAll.map(function(image){return image.multiply(10000).floor().toInt16()})
           .map(function(image){return image.clip(rectangle)});
print(lAll);

// There was a problem extracting pixel values. SampleRegions ignores pixels
// where at least one image has masked value. Thus, we need to unmask all the images
// and set footprint:false.
var timeSeries = ee.ImageCollection(lAll).toBands().unmask({value: -11000, sameFootprint: false});
// print(roi,'roi');
print(timeSeries, ' timeSeries');

// Map.addLayer(ee.Image(timeSeries.select(28)).clip(rectangle),{min:0,max:10000},'timeSeries');

// Get NDVI time series values as table
var sampledReg = ee.Image(timeSeries)
      .sampleRegions({
        // Get the sample from the polygons FeatureCollection.
        collection: pts4ts,
        properties: ['id','FrstTyp','Change'],
        // Set the scale to get Sentinel pixels in the polygons.
        scale: 30,
        tileScale: 8,
        geometries: true
      });

// Transform point coordinates into properties in the table.
var collection_with_latlon = sampledReg.map(function (feature) {
  var coordinates = feature.geometry().transform('epsg:4326').coordinates();
  return feature.set('lon', coordinates.get(0), 'lat', coordinates.get(1));
});

// Export table as csv
Export.table.toDrive({
    collection: collection_with_latlon,
    folder: folder1,
    description: 'verifiedpts4timeSeries_Disturb',
    fileFormat: 'CSV',
});

// Export first image just to check it
Export.image.toDrive({
                  image: ee.Image(timeSeries.select(1)).clip(rectangle),
                  description: '1eraImagenTimeSeriesBfastGEE',
                  scale: resolution,
                  folder: folder1,
                  crs: 'EPSG:4326',
                  maxPixels:1e12,
                  region: rectangle
                });