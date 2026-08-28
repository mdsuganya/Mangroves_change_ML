function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask).divide(10000);
}


var S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds (Bhitar);
  
//Add spectral indices
var addIndicesS2 = function(img){
  var ndvi = img.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var bndvi = img.normalizedDifference(['B8', 'B2']).rename('BNDVI');
  var ndwi = img.normalizedDifference(['B3', 'B8']).rename('NDWI');
  var mndwi = img.normalizedDifference(['B3', 'B11']).rename('MNDWI');
  var ndmi = img.normalizedDifference(['B8', 'B11']).rename('NDMI');
  var ndbi = img.normalizedDifference(['B11', 'B8']).rename('NDBI');
  var mdi1 = img.expression('(NIR-SWIR2)/(SWIR2)',{
    'NIR' :img.select('B8'),
    'SWIR2' :img.select('B12'),
  }).rename('MDI1');
   var mdi2 = img.expression('(NIR-SWIR1)/(SWIR1)',{
    'NIR' :img.select('B8'),
    'SWIR1' :img.select('B11'),
  }).rename('MDI2');
  var sipi = img.expression('(NIR-B)/(NIR-R)',{
    'NIR':img.select('B8'),
    'R' :img.select('B4'),
    'B' :img.select('B2')
  }).rename('SIPI');
 var savi = img.expression('1.25*(NIR-R)/(NIR+R+0.2)',{
    'NIR' :img.select('B8'),
    'R':img.select('B4')
  }).rename('SAVI');
  var sr = img.select('B8').divide(img.select('B4')).rename('SR');
  var ratio118 = img.select('B11').divide(img.select('B8')).rename('R118');
  var ratio411 = img.select('B4').divide(img.select('B11')).rename('R411');
  var gcvi = img.expression('(NIR/GREEN)-1',{
    'NIR':img.select('B8'),
    'GREEN':img.select('B3')
  }).rename('GCVI');
  var evi = img.expression('2.5*((NIR - R)/(NIR+6*R-7.5*B+1))',{
  'NIR':img.select('B8'),
    'R' :img.select('B4'),
    'B' :img.select('B2')
  }).rename('EVI');
  var arvi2 = img.expression('-0.18+1.17*(NIR-R)/(NIR+R)',{
    'NIR' :img.select('B8'),
    'R' :img.select('B4')
  }).rename('ARVI2');
  var bwdrvi = img.expression('(0.1*NIR-B)/(0.1*NIR+B)',{
     'NIR':img.select('B8'),
      'B' :img.select('B2')
  }).rename('BWDRVI');
  var ccci = img.expression('((NIR-RE)/(NIR+RE))/((NIR-R)/(NIR+R))',{
    'NIR':img.select('B8'),
      'RE' :img.select('B5'),
      'R'  :img.select('B4')
  }).rename('CCCI');
  var cig = img.expression('(NIR/G) - 1',{
    'NIR':img.select('B8'),
      'G' :img.select('B3')
  }).rename('CIG');
  var cire = img.expression('(NIR/G) - 1',{
    'NIR':img.select('B8'),
      'G' :img.select('B5')
  }).rename('CIRE');
  var gdvi = img.expression('(NIR-G)',{
    'NIR':img.select('B8'),
      'G' :img.select('B3')
  }).rename('GDVI');
  var evi2 = img.expression('2.4*(NIR-R)/(NIR+R+1)',{
    'NIR' :img.select('B8'),
    'R' :img.select('B4')
  }).rename('EVI2');
  var gvmi = img.expression('(NIR+0.1)-(SWIR+0.02)/(NIR+0.1)+(SWIR+0.02)',{
    'NIR' :img.select('B8'),
    'SWIR' :img.select('B12')
  }).rename('GVMI');
  var gari = img.expression('(NIR)-(G)-(B-R)/(NIR)-(G)+(B-R)',{
    'NIR' :img.select('B8'),
    'G' :img.select('B3'),
    'B' :img.select('B5'),
    'R' :img.select('B4')
  }).rename('GARI');
  var gli = img.expression('(2*G-R-B)/(2*G+R+B)',{
    'G' :img.select('B3'),
    'B' :img.select('B5'),
    'R' :img.select('B4')
  }).rename('GLI');
  var afri1 = img.expression('(VRE-0.66*SWIR)/(VRE+0.66*SWIR)',{
    'VRE' :img.select('B8A'),
    'SWIR' :img.select('B11')
  }).rename('AFRI1');
  var afri2 = img.expression('(VRE-0.5*SWIR)/(VRE+0.5*SWIR)',{
    'VRE' :img.select('B8A'),
    'SWIR' :img.select('B12')
  }).rename('AFRI2');
  var ari = img.expression('(1/G)-(1/VRE)',{
    'VRE' :img.select('B5'),
    'G' :img.select('B3')
  }).rename('ARI');
  var arvi = img.expression('NIR-(R-(B-R))/NIR+(R-(B-R))',{
    'NIR' :img.select('B8'),
    'R' :img.select('B4'),
    'B' :img.select('B2')
  }).rename('ARVI');
  var avi = img.expression('2*(VRE-R)',{
    'VRE' :img.select('B8A'),
    'R' :img.select('B4'),
  }).rename('AVI');
  var bri = img.expression('(1/G-1/VRE)/VRE1',{
    'G' :img.select('B3'),
    'VRE' :img.select('B5'),
    'VRE1' :img.select('B6')
  }).rename('BRI');
  var cri550 = img.expression('(1/B-1/G)',{
    'B' :img.select('B2'),
    'G' :img.select('B3'),
  }).rename('CRI550');
  var cvi = img.expression('(NIR*R)/(G**2)',{
    'NIR' :img.select('B8'),
    'R' :img.select('B4'),
    'G' :img.select('B3'),
  }).rename('CVI');
  var gndvi = img.expression('(NIR-G)/(NIR+G)',{
    'NIR' :img.select('B8'),
    'G' :img.select('B3')
  }).rename('GNDVI');
  var gndvi2 = img.expression('(VRE-G)/(VRE+G)',{
    'VRE' :img.select('B7'),
    'G' :img.select('B3')
  }).rename('GNDVI2');
  var ndii = img.expression('(NIR-SWIR)/(NIR+SWIR)',{
    'NIR' :img.select('B8'),
    'SWIR' :img.select('B11')
  }).rename('NDII');
  var osavi = img.expression('1.16*(NIR-R)/(NIR+R+0.16)',{
    'NIR' :img.select('B8'),
    'R' :img.select('B4')
  }).rename('OSAVI');
  var tcari = img.expression('3*((VRE-R)-0.2*(VRE-G)*(VRE/R))',{
    'VRE' :img.select('B5'),
    'R' :img.select('B4'),
    'G' :img.select('B3')
  }).rename('TCARI'); 

return img
  .addBands(ndvi)
  .addBands(bndvi)
  .addBands(ndwi)
  .addBands(mndwi)
  .addBands(ndmi)
  .addBands(ndbi)
  .addBands(mdi1)
  .addBands(mdi2)
  .addBands(sipi)
  .addBands(savi)
  .addBands(sr)
  .addBands(ratio118)
  .addBands(ratio411)
  .addBands(gcvi)
  .addBands(evi)
  .addBands(arvi2)
  .addBands(bwdrvi)
  .addBands(ccci)
  .addBands(cig)
  .addBands(cire)
  .addBands(gdvi)
  .addBands(evi2)
  .addBands(gvmi)
  .addBands(gari)
  .addBands(gli)
  .addBands(afri1)
  .addBands(afri2)
  .addBands(ari)
  .addBands(arvi)
  .addBands(avi)
  .addBands(bri)
  .addBands(cri550)
  .addBands(cvi)
  .addBands(gndvi)
  .addBands(gndvi2)
  .addBands(ndii)
  .addBands(osavi)
  .addBands(tcari)
};

  var filtered = S2
  .filterDate ('2024-01-01', '2024-12-30')
  .filter (ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .filterBounds (Bhitar)
  .map(maskS2clouds)
  .map(addIndicesS2);
  
var composite = filtered.median ();

var clipped = composite.clip(Bhitar);

var visbandsfcc = {bands: ['B6', 'B4', 'B3'], min : 0, max :0.3};
              
Map.addLayer(clipped, visbandsfcc, 'Sentinel 2024');
Map.centerObject (Bhitar, 10);
  print(composite);

// Map.addLayer(composite.clip(Bhitar),{bands:['B6', 'B4', 'B3'], min: 0, max: 0.3}, 'Sentinel 2024');

var srtmclip = SRTM.clip(Bhitar);
var elevationMask = srtmclip.lt(80);
var NDVIMask = composite.select('NDVI').gt(0.40);
var MNDWIMask = composite.select('MNDWI').gt(-0.80);
var NDBIMask = composite.select('NDBI').gt(-0.9);
var NDMIMask = composite.select('NDMI').gt(0.28);

//Apply the masks
var compositeNew = composite
                        .updateMask(NDVIMask)
                        .updateMask(MNDWIMask)
                        .updateMask(NDBIMask)
                        // .updateMask(NDMIMask)
                        .updateMask(elevationMask)
var compositefloat = composite.toFloat()
        
                      
//2.7) Display results
//////////////////////                       
//Select bands and parameters for visualization
var visPar = {bands:['B6','B4','B3'], min: 0, max: 0.3}; 

//Add layer to map
Map.addLayer(compositeNew.clip(Bhitar), visPar, 'Sentinel Composite 2024');

//construct Random Model

var classes = Mangroves.merge(nonmangrove).merge(River)

var bands = ['B12', 'B11', 'B8', 'B7', 'B6', 'B5', 'B4', 'B3', 'B2', 'NDVI','NDMI','NDBI', 'MNDWI', 'SR','GCVI']

var image = composite.select(bands).clip(Bhitar)

var samples = image.sampleRegions({
    collection: classes,
    properties: ['landcover'],
    scale: 30,
    tileScale: 16,
    }).randomColumn('random');
    
var split = 0.8;
var training = samples.filter(ee.Filter.lt('random', split));
var testing = samples.filter(ee.Filter.gte('random', split));

    print('Samples n =', samples.aggregate_count('.all'));
    print('Training n =', training.aggregate_count('.all'));
    print('Testing n =', testing.aggregate_count('.all'));

//3.2) Begin Random Forest Classification

//.smileRandomForest is used to run the model. Here we run the model using 100 trees
// and 5 randomly selected predictors per split ("(100,5)")

    var classifier = ee.Classifier.smileRandomForest(1000,5).train({
    features: training.select('B12', 'B11', 'B8', 'B7', 'B6', 'B5', 'B4', 'B3', 'B2', 'NDVI','NDMI','NDBI', 'MNDWI', 'SR','GCVI', 'landcover'),
    classProperty: 'landcover',
    inputProperties: bands
    });
    
//3.3) Test the accuracy of the model

    var validation = testing.classify(classifier);
    var testAccuracy = validation.errorMatrix('landcover', 'classification');
    print('Validation error matrix RF: ', testAccuracy);
    print('Validation overall accuracy RF: ', testAccuracy.accuracy());

//3.4) Classify the Sentinel-2 composite using the Random Forest model
/////////////////////////////////////////////////////////////////

    var classifiedrf = image.select(bands)
                      .classify(classifier);
                      
//The model results may be "noisy". To reduce noise, create a mask to mask
// unconnected pixels

    var pixelcount = classifiedrf.connectedPixelCount(100, false);
    var countmask = pixelcount.select(0).gt(25)
//Mask the results to only display mangrove extent

    var classMask = classifiedrf.select('classification').gt(0)
    var classed = classifiedrf.updateMask(countmask).updateMask(classMask)

//3.5) Map results
////////////////////

//Add classification to map

Map.addLayer (classed, {min: 1, max: 3, palette:['green','blue','red']}, 'Mangrove Extent 2024');
var mangclip = classed.clip(Bhitar);
Map.addLayer (mangclip,{min: 1, max: 3, palette:['green','blue','red']}, 'Mangrove clip 2024');
                      
//2024 Mangrove Extent
//------------------
Export.image.toDrive({
  image: classed,
  description: '2024BhitarkanikaMangroveExtent',
  region: Bhitar,
  scale: 30,
  crs: 'EPSG:32645',
  maxPixels: 1e13
  }); 
