var geometry = 
    /* color: #d6c4c3 */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[4.645213482087542, 43.61416821722853],
          [4.645213482087542, 43.60913454258402],
          [4.652852413361956, 43.60913454258402],
          [4.652852413361956, 43.61416821722853]]], null, false);
//GEO441 Camargue Indices calculation
//ROI
/*
var geometry = 
    ee.Geometry.Polygon(
        [[[3.975781835937493, 43.93432837314937],
          [3.975781835937493, 43.08577150036554],
          [5.447949804687493, 43.08577150036554],
          [5.447949804687493, 43.93432837314937]]], null, false);
*/

//Pre setting for filter and image selection
var cloud_cover = 80    // to filter images above the threshold

var i = 20  // select the i-th image of the collection for mapping

//time
var startyear = 2015; 
var endyear = 2019; 
var startmonth = 1;
var endmonth = 12;
var startday = 1; //
var endday = 365; 

// SPACE

// AOI
var aoi_name = 'Camargue'; 

var aoi = geometry; 
print(aoi)

//Data
//Landsat SR collections

var oliCol = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR');
var etmCol = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR');
var tmCol = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR');

/*
var dataset = ee.ImageCollection('NASA/MEASURES/GFCC/TC/v3')
                  .filter(ee.Filter.date('2015-01-01', '2015-12-31'));
var treeCanopyCover = dataset.select('tree_canopy_cover');
var forestMask = treeCanopyCover.mean().gt(50).unmask().clip(aoi)
*/

// Define a collection filter
var colFilter = ee.Filter.and(
    ee.Filter.bounds(aoi), 
    // ee.Filter.calendarRange(startday, endday, 'day_of_year'),
    ee.Filter.calendarRange(startyear, endyear, 'year'),
    ee.Filter.calendarRange(startmonth, endmonth, 'month'),
    ee.Filter.lt('CLOUD_COVER', cloud_cover), 
    ee.Filter.lt('GEOMETRIC_RMSE_MODEL', 10),
    ee.Filter.or(
        ee.Filter.eq('IMAGE_QUALITY', 9),
        ee.Filter.eq('IMAGE_QUALITY_OLI', 9)));  
  
//Visualization

var visParams = {
  opacity: 1,
  bands: ['slope'],
  min: -50,
  max: 50,
  palette:
    ['8c510a', 'd8b365', 'f6e8c3', 'f5f5f5', 'd9f0d3', '7fbf7b', '1b7837']
};

var assets = require('users/gena/packages:palettes')

var VisParams_NDVI = {   min: 0.2, max: 0.8, 
                    palette: assets.colorbrewer.RdYlGn[9]}

var VisParams_NDMI = {   min: 0.2, max: 0.8, 
                    palette: assets.colorbrewer.RdYlGn[9]}

var VisParams_NDWI = {   min: 0.2, max: 0.8, 
                    palette: assets.colorbrewer.RdYlGn[9]}

var VisParams_BSI = {   min: -1, max: 1, 
                    palette: assets.colorbrewer.RdYlGn[9]}

var VisParams_EVI = {   min: 1, max: 4, 
                    palette: assets.colorbrewer.RdYlGn[9]}

var VisParams_IND1 = {   min: -1, max: 1, 
                    palette: assets.colorbrewer.RdYlGn[9]}

var VisParams_IND2 = {   min: -1, max: 1, 
                    palette: assets.colorbrewer.RdYlGn[9]}

Map.setOptions('SATELLITE');
Map.addLayer(aoi, {color: 'white'}, 'AOI', true);
Map.centerObject(aoi,12);


// PRE-PROCESSING 

// Harmonization of ETM+ spectral space to OLI spectral space 
// Band-respective coefficients are defined in the following dictionary with slope (slopes) and intercept (itcps) image constants. 
// Note that the y-intercept values are multiplied by 10,000 to match the scaling of USGS Landsat surface reflectance data.

var coefficients = {
  itcps: ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172])
             .multiply(10000),
  slopes: ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071])
};

// ETM+ and OLI band names for the same spectral response window are not equal and need to be standardized. 
// The following functions select only reflectance bands and the pixel_qa band from each dataset, 
// and renames them according to the wavelength range they represent.

// Function to get and rename bands of interest from OLI.
function renameOli(img) {
  return img.select(
      ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'pixel_qa'],
      ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'pixel_qa']);
}

// Function to get and rename bands of interest from ETM+.
function renameEtm(img) {
  return img.select(
      ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'pixel_qa'],
      ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'pixel_qa']);
}

// Finally, define the transformation function, which applies the linear model to ETM+ data, 
// casts data type as Int16 for consistency with OLI, 
// and reattaches the pixel_qa band for later use in cloud and shadow masking.

function etmToOli(img) {
  return img.select(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])
      .multiply(coefficients.slopes)
      .add(coefficients.itcps)
      .round()
      .toShort()
      .addBands(img.select('pixel_qa'));
}

// Cloud and shadow masking

// The following function uses the CFmask (Zhu et al., 2015) pixel_qa band 
// included with each Landsat USGS surface reflectance image 
// to set pixels identified as cloud and cloud shadow to null.

function fmask(img) {
  var cloudShadowBitMask = 1 << 3;
  var cloudsBitMask = 1 << 5;
  var qa = img.select('pixel_qa');
  var mask = qa.bitwiseAnd(cloudShadowBitMask)
                 .eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return img.updateMask(mask);
}


// Spectral index calculation

// The normalized burn ratio (NBR) spectral index represents the spectral history of a forested pixel affected by wildfire.
function calcNDWI(img) {
  return img.addBands(img.normalizedDifference(['Green', 'NIR']).rename('NDWI'));
}

// The normalized difference moisture index (NDMI) describes the vegetations's water stress level
function calcNDMI(img) {
  return img.addBands(img.normalizedDifference(['NIR', 'SWIR1']).rename('NDMI'));
}

// The normalized vegetation index (NDVI) estimates the density of vegetation on an area of land.
function calcNDVI(img) {
  return img.addBands(img.normalizedDifference(['NIR', 'Red']).rename('NDVI'));
}  

// The normalized bare soil index (BSI)
function calcBSI(img){
  var bsi = img.expression(
    '((Red+SWIR1)-(NIR+Blue))/((Red+SWIR1)+(NIR+Blue))', {
      'SWIR1': img.select('SWIR1'),
      'NIR': img.select('NIR'),
      'Red': img.select('Red'),
      'Blue': img.select('Blue')
})
  
  return img.addBands(bsi.rename('BSI'));
}

// The enhanced vegetation index EVI
function calcEVI(img) {
  var evi = img.expression(
    '2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)',
    {
        'red': img.select('Red'),    
        'nir': img.select('NIR'),    
        'blue': img.select('Blue')    
    });
  
  return img.addBands(evi.rename('EVI'));
}

// The new index LSWI(NDMI)-EVI+0.05
function calcIND1(img) {
  var IND1 = img.expression(
    '(nir - swir1) / (nir + swir1) - 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1) + 0.05',
    {
      'swir1': img.select('SWIR1'),
      'nir': img.select('NIR'),
      'red': img.select('Red'),
      'blue': img.select('Blue')
    });
  
  return img.addBands(IND1.rename('IND1'));
}

// The new index LSWI(NDMI)-NDVI+0.05
function calcIND2(img) {
  var IND2 = img.expression(
    '(nir - swir1) / (nir + swir1) - (nir - red) / (nir +  red) + 0.05',
    {
      'swir1': img.select('SWIR1'),
      'nir': img.select('NIR'),
      'red': img.select('Red'),
    });
  
  return img.addBands(IND2.rename('IND2'));
}



// PROCESSING 
// Define function to prepare OLI images.
function prepOli(img) {
  var orig = img;
  img = renameOli(img);
  img = fmask(img);
  img = calcNDWI(img);
  img = calcNDMI(img);
  img = calcNDVI(img);
  img = calcBSI(img);
  img = calcEVI(img);
  img = calcIND1(img);
  img = calcIND2(img);
  return ee.Image(img.copyProperties(orig, orig.propertyNames()));
}

// Define function to prepare ETM+ images.
function prepEtm(img) {
  var orig = img;
  img = renameEtm(img);
  img = fmask(img);
  img = calcNDWI(img);
  img = calcNDMI(img);
  img = calcNDVI(img);
  img = calcBSI(img);
  img = calcEVI(img);
  img = calcIND1(img);
  img = calcIND2(img);
  return ee.Image(img.copyProperties(orig, orig.propertyNames()));
}

///////////////////////////////////////////////////////
// Combine collections

// Filter the collections and map the prepImg function over all images. 
oliCol = oliCol.filter(colFilter).map(prepOli);
etmCol = etmCol.filter(colFilter).map(prepEtm);
tmCol = tmCol.filter(colFilter).map(prepEtm);

// Merge the collections.
var col = oliCol.merge(etmCol).merge(tmCol);
// print(col.first())



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MAPPING 

// stick a selected image from the collection to the map

print('observation date', ee.Date(ee.Image(col.toList(1,i).get(0)).get('system:time_start')))

Map.addLayer(ee.Image(col.toList(1,i).get(0)).select('NDWI'), VisParams_NDWI, 'NDWI', true)
Map.addLayer(ee.Image(col.toList(1,i).get(0)).select('NDMI'), VisParams_NDMI, 'NDMI', true)
Map.addLayer(ee.Image(col.toList(1,i).get(0)).select('NDVI'), VisParams_NDVI, 'NDVI', true)
Map.addLayer(ee.Image(col.toList(1,i).get(0)).select('BSI'), VisParams_BSI, 'BSI', true)
Map.addLayer(ee.Image(col.toList(1,i).get(0)).select('EVI'), VisParams_EVI, 'EVI', true)
Map.addLayer(ee.Image(col.toList(1,i).get(0)).select('IND1'), VisParams_IND1, 'IND1', true)
Map.addLayer(ee.Image(col.toList(1,i).get(0)).select('IND2'), VisParams_IND2, 'IND2', true)


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ANALYSIS 

// Map a region reduce function over the image collection that calculates the median of all pixels intersecting the AOI 
// and sets the result as an image property.

// NDWI
var obsNDWI = col.map(function(img) {
  var obs = img.select('NDWI').reduceRegion(
      {geometry: aoi, reducer: ee.Reducer.median(), scale: 30});
  return img.set('NDWI_prop', obs.get('NDWI'));
});
print(obsNDWI)

// Plotting all observations as a scatter plot where sensor is discerned by color
var chartObsNDWI =
    ui.Chart.feature.groups(obsNDWI, 'system:time_start', 'NDWI_prop', 'SATELLITE')
        .setChartType('ScatterChart')
        .setSeriesNames(['ETM+', 'OLI']) // 'TM', 
        .setOptions({
          title: 'NDWI Observations',
          colors: ['f8766d', '00ba38', '619cff'],
          hAxis: {title: 'Date'},
          vAxis: {title: 'NDWI'},
          pointSize: 6,
          dataOpacity: 0.5
        });
print(chartObsNDWI);


// NDMI
var obsNDMI = col.map(function(img) {
  var obs = img.select('NDMI').reduceRegion(
      {geometry: aoi, reducer: ee.Reducer.median(), scale: 30});
  return img.set('NDMI_prop', obs.get('NDMI'));
});

// print(obsNDMI.first())

// Plotting all observations as a scatter plot where sensor is discerned by color
var chartObsNDMI =
    ui.Chart.feature.groups(obsNDMI, 'system:time_start', 'NDMI_prop', 'SATELLITE')
        .setChartType('ScatterChart')
        .setSeriesNames(['ETM+', 'OLI']) // 'TM', 
        .setOptions({
          title: 'NDMI Observations',
          colors: ['f8766d', '00ba38', '619cff'],
          hAxis: {title: 'Date'},
          vAxis: {title: 'NDMI'},
          pointSize: 6,
          dataOpacity: 0.5
        });
print(chartObsNDMI);


// NDVI
var obsNDVI = col.map(function(img) {
  var obs = img.select('NDVI').reduceRegion(
      {geometry: aoi, reducer: ee.Reducer.median(), scale: 30});
  return img.set('NDVI_prop', obs.get('NDVI'));
});

// Plotting all observations as a scatter plot where sensor is discerned by color
var chartObsNDVI =
    ui.Chart.feature.groups(obsNDVI, 'system:time_start', 'NDVI_prop', 'SATELLITE')
        .setChartType('ScatterChart')
        .setSeriesNames(['ETM+', 'OLI']) // 'TM', 
        .setOptions({
          title: 'NDVI Observations',
          colors: ['f8766d', '00ba38', '619cff'],
          hAxis: {title: 'Date'},
          vAxis: {title: 'NDVI'},
          pointSize: 6,
          dataOpacity: 0.5
        });
print(chartObsNDVI);


// BSI
var obsBSI = col.map(function(img) {
  var obs = img.select('BSI').reduceRegion(
      {geometry: aoi, reducer: ee.Reducer.median(), scale: 30});
  return img.set('BSI_prop', obs.get('BSI'));
});

// Plotting all observations as a scatter plot where sensor is discerned by color
var chartObsBSI =
    ui.Chart.feature.groups(obsBSI, 'system:time_start', 'BSI_prop', 'SATELLITE')
        .setChartType('ScatterChart')
        .setSeriesNames(['ETM+', 'OLI'])    // 'TM', 
        .setOptions({
          title: 'BSI Observations',
          colors: ['f8766d', '00ba38', '619cff'],
          hAxis: {title: 'Date'},
          vAxis: {title: 'BSI'},
          pointSize: 6,
          dataOpacity: 0.5
        });
print(chartObsBSI);

var chartMedianComp = ui.Chart.image
                          .series({
                            imageCollection: col,
                            region: aoi,
                            reducer: ee.Reducer.median(),
                            scale: 30,
                            xProperty: 'system:time_start',
                          })
                          .setSeriesNames(['NDWI Median'])
                          .setOptions({
                            title: 'Intra-annual Median',
                            colors: ['619cff'],
                            hAxis: {title: 'Date'},
                            vAxis: {title: 'NDWI'},
                            lineWidth: 6
                          });
print(chartMedianComp);

//EVI
var obsEVI = col.map(function(img) {
  var obs = img.select('EVI').reduceRegion(
      {geometry: aoi, reducer: ee.Reducer.median(), scale: 30});
  return img.set('EVI_prop', obs.get('EVI'));
});

// Plotting all observations as a scatter plot where sensor is discerned by color
var chartObsEVI =
    ui.Chart.feature.groups(obsEVI, 'system:time_start', 'EVI_prop', 'SATELLITE')
        .setChartType('ScatterChart')
        .setSeriesNames(['ETM+', 'OLI'])    // 'TM', 
        .setOptions({
          title: 'EVI Observations',
          colors: ['f8766d', '00ba38', '619cff'],
          hAxis: {title: 'Date'},
          vAxis: {title: 'EVI'},
          pointSize: 6,
          dataOpacity: 0.5
        });
print(chartObsEVI);

// IND1
var obsIND1 = col.map(function(img) {
  var obs = img.select('IND1').reduceRegion(
      {geometry: aoi, reducer: ee.Reducer.median(), scale: 30});
  return img.set('IND1_prop', obs.get('IND1'));
});

// Plotting all observations as a scatter plot where sensor is discerned by color
var chartObsIND1 =
    ui.Chart.feature.groups(obsIND1, 'system:time_start', 'IND1_prop', 'SATELLITE')
        .setChartType('ScatterChart')
        .setSeriesNames(['ETM+', 'OLI'])    // 'TM', 
        .setOptions({
          title: 'IND1 Observations',
          colors: ['f8766d', '00ba38', '619cff'],
          hAxis: {title: 'Date'},
          vAxis: {title: 'IND1'},
          pointSize: 6,
          dataOpacity: 0.5
        });
print(chartObsIND1);

// IND2
var obsIND2 = col.map(function(img) {
  var obs = img.select('IND2').reduceRegion(
      {geometry: aoi, reducer: ee.Reducer.median(), scale: 30});
  return img.set('IND2_prop', obs.get('IND2'));
});

// Plotting all observations as a scatter plot where sensor is discerned by color
var chartObsIND2 =
    ui.Chart.feature.groups(obsIND2, 'system:time_start', 'IND2_prop', 'SATELLITE')
        .setChartType('ScatterChart')
        .setSeriesNames(['ETM+', 'OLI'])    // 'TM', 
        .setOptions({
          title: 'IND2 Observations',
          colors: ['f8766d', '00ba38', '619cff'],
          hAxis: {title: 'Date'},
          vAxis: {title: 'IND2'},
          pointSize: 6,
          dataOpacity: 0.5
        });
print(chartObsIND2);
