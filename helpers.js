
//===============================================================================================================================================
// auxiliars functions to get images collections and vegetation index
// getSr - conjunto de funções para buscar imagens de acordo com o sensor
// indexSr - conjunto de funções usadas para gerar índices espectrais
// hantModel - modelo harmonico para suavizar séries temporais
// getArea - conjunto de funções para realizar cálculos de área


exports.getSr = {
  
    filterCollectionS1: function (instrumentMode, roi, asc_desc) {
        
        var col = ee.ImageCollection('COPERNICUS/S1_GRD')
            //.filter(ee.Filter.listContains('transmitterReceiverPolarisation', transRecivPol))
            .filter(ee.Filter.eq('instrumentMode', instrumentMode))
            .filterBounds(roi)
            .map(function (image) {
                image = image.clip(
                    image.geometry().buffer(-50).simplify(1)
                );
                return image;
            });
        
        if (asc_desc == 'asc') {
            col = col.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));
        } else {
            col = col.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));
        }
        
        return col;
        
    },
    
    filterCollectionS2: function (roi, t0, t1, cloud) {
      
        function removerNuvens (image){
          
          var cloud = ee.Image(image.select(['aerosols']));
          return image.mask(cloud.lt(2000));
          
        }

        var names = [
            'B1','B2', 'B3', 'B4', 'B5', 'B6', 'B7','B8', 'B8A','B11', 'B12'
        ];
        var newNames = [
            'aerosols','blue', 'green', 'red', 'redEdge1', 'redEdge2', 
            'redEdge3', 'nir', 'redEdge4', 'swir1', 'swir2'
        ];
        
        var col = ee.ImageCollection('COPERNICUS/S2')
            .filterBounds(roi)
            .filterDate(t0, t1)
            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloud))
            .map(function (image) { return image.clip(roi) })
            .select(names, newNames)
            .map(removerNuvens);
        
        return col;
    },
    
    filterCollectionL8: function (roi, t0, t1, cloud) {
      
        function removeCloud(image) {
          // Bits 3 and 5 are cloud shadow and cloud, respectively.
          var cloudShadowBitMask = (1 << 3);
          var cloudsBitMask = (1 << 5);
          // Get the pixel QA band.
          var qa = image.select('pixel_qa');
          // Both flags should be set to zero, indicating clear conditions.
          var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                         .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
          return image.updateMask(mask);
        }
        
        var names = [
            'B1','B2', 'B3', 'B4', 'B5', 'B6', 'B7','B10', 'B11','pixel_qa'
        ];
        var newNames = [
            'ultra_blue','blue', 'green', 'red', 'nir', 'swir1', 
            'swir2', 'brightness_temp1', 'brightness_temp2', 'pixel_qa'
        ];
        
        var col = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
            .filterBounds(roi)
            .filterDate(t0, t1)
            .filter(ee.Filter.lt('CLOUD_COVER', cloud))
            .map(function (image) { return image.clip(roi) })
            .select(names, newNames)
            .map(removeCloud);
        
        return col;      
    },

    filterCollectionL5: function (roi, t0, t1, cloud) {
      
        function removeCloud (image){
          
          var qa = image.select('pixel_qa');
          var cloud = qa.bitwiseAnd(1 << 5)
                    .and(qa.bitwiseAnd(1 << 7))
                    .or(qa.bitwiseAnd(1 << 3));
              
          var mask2 = image.mask().reduce(ee.Reducer.min());
          return image.updateMask(cloud.not()).updateMask(mask2);
          
        }
        
        var bands = {
          names: ['B1','B2', 'B3', 'B4', 'B5', 'B6','B7', 'sr_atmos_opacity', 'sr_cloud_qa'],
          newNames: ['blue', 'green', 'red','nir', 'swir1','temp','swir2','sr_atmos_opacity', 'sr_cloud_qa' ]
        };
        
        var col = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
            .filterBounds(roi)
            .filterDate(t0, t1)
            .filter(ee.Filter.lt('CLOUD_COVER', cloud))
            .map( function(i) { return i.clip(roi) })
            .map(removeCloud)
            .select(bands.names, bands.newNames);
        
        return col;      
    }
};

exports.miningIndex = {
  
    // Bare Soil Index (BSI) is a numerical indicator that combines blue, red, near infrared and short wave 
    // infrared spectral bands to capture soil variations. These spectral bands are used in a normalized manner. 
    // The short wave infrared and the red spectral bands are used to quantify the soil mineral composition, while 
    // the blue and the near infrared spectral bands are used to enhance the presence of vegetation (GU, 2019).
    
    addBSI: function(img){

        var BSI = img.expression('(red + swir1)-(nir + blue) / (red + swir1) + (nir + blue)', {
            'nir': img.select('nir'),
            'swir1': img.select('swir1'),
            'blue': img.select('blue'),
            'red': img.select('red')
        });
        
        return img.addBands(BSI.rename('BSI'));
    },
    
    // Iron Oxide
    addIo: function(img){
        var index = img.expression('red / blue', {
            'red': img.select('red'),
            'blue': img.select('blue')
        });
        return img.addBands(index.rename('io'));
    },
    
    // Fe3+ - B04 / B03;
    addFe3: function(img){
        var index = img.expression('red / green', {
            'red': img.select('red'),
            'green': img.select('green')
        });
        return img.addBands(index.rename('fe3'));
    },
    
    //Fe2+ - B12 / B08 + B03 / B04;
    addFe2: function(img){
        var index = img.expression('(swir2 / nir) + (green / red)', {
            'swir2': img.select('swir2'),
            'nir': img.select('nir'),
            'green': img.select('green'),
            'red': img.select('red')
        });
        return img.addBands(index.rename('fe2'));
    },

    // Simple Ratio SWIRI/NIR Ferrous Minerals (abbrv. SRSWIRI/NIR)
    addFerrousMinerals: function(img){
        var index = img.expression('(swir1 / nir)', {
            'swir1': img.select('swir1'),
            'nir': img.select('nir')
        });
        return img.addBands(index.rename('ferrousMinerals'));
    },
    
};

exports.indexSr = {

    addSavi: function (image) {
    
        var savi = image.expression('(1 + L) * float(nir - red)/(nir + red + L)', {
            'nir': image.select('nir'),
            'red': image.select('red'),
            'L': 0.9
        });
    
        return image.addBands(savi.rename("savi"));
    
    },
    
    addEviL: function (image) {
    
        var evi = image.expression('2.5 * ( (nir-red)-(nir+2.4*red+1) )', {
            'nir': image.select('nir'),
            'red': image.select('red'),
        });
    
        return image.addBands(evi.rename("evi"));
    
    },

    addEviSe: function (image) {
    
        var evi = image.expression('2.5 * ( (nir-red)-(nir+2.4*red+1) )', {
            'nir': image.select('nir'),
            'red': image.select('red'),
        });
    
        return image.addBands(evi.rename("evi"));
    
    },
    
    addNdvi: function (image){
        var ndvi = image.normalizedDifference(['nir','red']);
        return image.addBands(ndvi.rename('ndvi'));
    
    },
    
    addNdwi: function (image) {
    
        var ndwi = image.normalizedDifference(['nir','swir1']);
        return image.addBands(ndwi.rename('ndwi'));
    
    },
    
    addSti: function (image) {
        var sti = image.select('swir1').divide(image.select('swir2'));
        return image.addBands(sti.rename('sti'));
    },
    
    addNdti: function (image) {
        var ndti = image.normalizedDifference(['swir1','swir2']);
        return image.addBands(ndti.rename('ndti'));      
    },

    addLac: function (image) {
        var lac = image.expression('blue-redEdge1', {
            'blue':image.select('blue'),
            'redEdge1':image.select('redEdge1')
        });
        return image.addBands(lac.rename('lac'));      
    },
    
    
    
    addMetrics: function (collection, image, band) {
      
      var max = collection.select(band).max().rename(band + '_max');
      var min = collection.select(band).min().rename(band + '_min');
      
      return image.addBands([max, min]);
      
    }

};

exports.hantModel = {
  
  //NDVIt = β0 + β1t + β2cos(2πωt) + β3sin(2πωt) + …
  
    getfitValues: function (n_cycles, dependent, independents, timeDif, collection) { // independents = [t, constant]
      
        var harmonicTrend = this.getCof(n_cycles, dependent, independents, timeDif, collection);
        
        var harmonicTrendCoefficients = harmonicTrend[1].select('coefficients')
          .arrayProject([0])
          .arrayFlatten([ harmonicTrend[2] ]);
          
        function fitValues(image){
          return image.addBands(image.select(harmonicTrend[2])
                .multiply(harmonicTrendCoefficients)
                .reduce('sum')
                .rename('fitted'));
        }
        
        var fittedHarmonic = harmonicTrend[0].map(fitValues);
        
        return [fittedHarmonic, harmonicTrendCoefficients];
          
    },
    
    getCof: function (n_cycles, dependent, independents, timeDif, collection) {
        
        var harmonicFrequencies = ee.List.sequence(1, n_cycles);
        
        function bandNames(base, list) { 
            return ee.List(list).map( function(i) {
                return ee.String(base).cat(ee.Number(i).int());
          });
        }
        
        var cosName = bandNames('cos_', harmonicFrequencies);
        var sinName = bandNames('sin_', harmonicFrequencies);
        
        independents = cosName.cat(sinName).cat(independents);
        
        function addDependents (image) {
            
            var years = image.date().difference(timeDif, 'year');
            var timeRadians = ee.Image(years.multiply(2 * Math.PI)).rename('t');
            var constant = ee.Image(1);
            
            return image.addBands(constant).addBands(timeRadians.float());
            
        }
        
        function addHarmonics (freqs) {
            return function (image) {
              
                var frequence = ee.Image.constant(freqs);
                var time = ee.Image(image).select('t');
                
                var sin = time.multiply(frequence).sin().rename(sinName);
                var cos = time.multiply(frequence).cos().rename(cosName);
                
                return image.addBands(cos).addBands(sin);
            };
        }
        
        var harmonicCollection = collection.map(addDependents).map( addHarmonics(harmonicFrequencies) );
        
        var harmonicTrend = harmonicCollection
          .select(independents.add(dependent))
          .reduce(ee.Reducer.linearRegression(independents.length(), 1));
        
        return [harmonicCollection, harmonicTrend, independents];
    }    
};

exports.getArea = {
  
  image:null,
  boundaries:null,
  prop:null,
  scale:null,
  
  loadData: function (boundaries, prop, image, scale) {
    this.image = image;
    this.boundaries = boundaries;
    this.prop = prop;
    this.scale = scale;
  },
  
  calcular: function (image, scale, boundaries) {
    
    var pxArea = ee.Image.pixelArea().divide(10000).multiply(image);
    
    var calc = pxArea.reduceRegion({
        reducer: ee.Reducer.sum(),
        scale:scale,
        geometry:boundaries,
        maxPixels: 1e13
    });
    
    return calc;
  },
 
  limitIdList: function (boundaries, prop) {
    
    var bounds = boundaries.toList( boundaries.size() );
    var idList = bounds.map( function (i) {
        return ee.Feature(i).get(prop);
    });
    
    return idList;
  },
  
  createListArea: function (boundaries, prop, image, scale) {
    
    function calcular (image, scale, boundaries) {
      
      var pxArea = ee.Image.pixelArea().divide(1000000).multiply(image);
      
      var calc = pxArea.reduceRegion({
          reducer: ee.Reducer.sum(),
          scale:scale,
          geometry:boundaries,
          maxPixels: 1e13
      });
      
      return calc;      
    }
    
    var id = this.limitIdList(boundaries, prop);
    
    var list = id.map( function (element) {
       
        var elementId = boundaries.filterMetadata(prop, 'equals', element);
        var img = image.clip( elementId ); 
        
        var area = calcular( img, scale, elementId.geometry() );
          
        return ee.Feature(null, area).set(prop, element);
        
    } );
      
    return list;
    
  },
  
  init: function (boundaries, prop, image, scale) {
    
    this.loadData(boundaries, prop, image, scale);
    return this.createListArea(this.boundaries, this.prop, this.image, this.scale);
    
  }
  
};
  
