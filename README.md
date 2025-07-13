ee.Authenticate()
ee.Initialize(project='my-project-dissertation-464420')

import ee
import folium # For interactive mapping in a Jupyter/Colab environment
import ipywidgets as widgets
from ipyleaflet import Map, TileLayer, WidgetControl, LayersControl

# Initialize Earth Engine (assuming it's already initialized in a previous cell with a project)
# The previous cell successfully initialized Earth Engine, so we don't need to do it again here.
# try:
#     ee.Initialize()
#     print("Earth Engine initialized successfully!")
# except ee.EEException as e:
#     print(f"Earth Engine initialization failed: {e}")
#     print("Please authenticate your GEE account. If you're running this for the first time,")
#     print("you might need to run 'earthengine authenticate' in your terminal.")
    # Exit or handle the error appropriately if initialization fails
    # exit()

# Define Rimba Raya AOI as a polygon
aoi = ee.Geometry.Polygon([
  [
    [112.020, -2.5],
    [112.470, -2.5],
    [112.470, -3.35],
    [112.020, -3.35],
    [112.020, -2.5]  # closing the ring
  ]
])

# --- Functions for Landsat Preprocessing ---

# Function to apply scale factors to Landsat Collection 2, Level 2 images (SR and ST bands).
def apply_scale_factors(image):
  optical_bands = image.select('SR_B.*').multiply(0.0000275).add(-0.2)
  thermal_bands = image.select('ST_B.*').multiply(0.00341802).add(149.0)
  return image.addBands(optical_bands, None, True).addBands(thermal_bands, None, True)

# Function to mask clouds and shadows in Landsat 8 Collection 2, Level 2 images (QA_PIXEL band).
def cloud_mask_l8(image):
  qa_pixel = image.select('QA_PIXEL')
  # Bits 3 and 5 are cloud shadow and cloud, respectively.
  cloud_shadow_bitmask = 1 << 3
  clouds_bitmask = 1 << 5
  # Both flags should be zero, indicating clear conditions.
  mask = qa_pixel.bitwiseAnd(cloud_shadow_bitmask).eq(0).And(
         qa_pixel.bitwiseAnd(clouds_bitmask).eq(0))
  return image.updateMask(mask)

# Function to mask clouds and shadows in Landsat 5 Collection 2, Level 2 images (QA_PIXEL band).
def cloud_mask_l5(image):
  qa_pixel = image.select('QA_PIXEL')
  # Bits 3 and 5 are cloud shadow and cloud, respectively.
  cloud_shadow_bitmask = 1 << 3
  clouds_bitmask = 1 << 5
  # Both flags should be zero, indicating clear conditions.
  mask = qa_pixel.bitwiseAnd(cloud_shadow_bitmask).eq(0).And(
         qa_pixel.bitwiseAnd(clouds_bitmask).eq(0))
  return image.updateMask(mask)

# Function to calculate NDVI for Landsat imagery.
# Assumes bands have been renamed to 'NIR' and 'Red'.
def calculate_ndvi(image):
  ndvi = image.normalizedDifference(['NIR', 'Red']).rename('NDVI')
  return image.addBands(ndvi)

# --- Define Time Period ---
start_year = 2008
end_year = 2025 # Data up to June 30, 2025

# Define yearsList (ee.List)
years_list = ee.List.sequence(start_year, end_year)
# print('Raw yearsList (EE Object):', years_list.getInfo()) # DEBUG with .getInfo()

# --- Load and Preprocess Landsat Data for the entire period ---
all_landsat = ee.ImageCollection([]) # Empty collection to start

# Landsat 5 (Operational 1984-03-01 to 2012-05-05 for C2L2)
collection_l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2') \
  .filterDate(str(start_year) + '-01-01', '2012-12-31') \
  .filterBounds(aoi)

# print('L5 Raw collection size:', collection_l5.size().getInfo()) # DEBUG
collection_l5 = collection_l5 \
  .map(apply_scale_factors) \
  .map(cloud_mask_l5) \
  .select(['SR_B4', 'SR_B3', 'QA_PIXEL'], ['NIR', 'Red', 'QA_PIXEL'])
# print('L5 After processing and band selection size:', collection_l5.size().getInfo()) # DEBUG
all_landsat = all_landsat.merge(collection_l5)

# Landsat 8 (Operational 2013-04-11 to present for C2L2)
collection_l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
  .filterDate('2013-01-01', str(end_year) + '-06-30') \
  .filterBounds(aoi)

# print('L8 Raw collection size:', collection_l8.size().getInfo()) # DEBUG
collection_l8 = collection_l8 \
  .map(apply_scale_factors) \
  .map(cloud_mask_l8) \
  .select(['SR_B5', 'SR_B4', 'QA_PIXEL'], ['NIR', 'Red', 'QA_PIXEL'])
# print('L8 After processing and band selection size:', collection_l8.size().getInfo()) # DEBUG
all_landsat = all_landsat.merge(collection_l8)

# print('Merged Landsat collection size (before NDVI calculation):', all_landsat.size().getInfo()) # DEBUG

# Calculate NDVI for the merged collection and select only the NDVI band
ndvi_collection = all_landsat.map(calculate_ndvi).select('NDVI')
# print('NDVI collection size (after NDVI calculation):', ndvi_collection.size().getInfo()) # DEBUG


# --- Visualization Parameters ---

# Visualization for individual year's NDVI (if you want to see them)
ndvi_vis = {
  'min': 0.1,   # Typical min for dense vegetation
  'max': 0.8,   # Typical max for dense vegetation
  'palette': [
    'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163',
    '99B718', '74A901', '66A000', '529400', '3e8601'
  ]
}

# Visualization for NDVI difference (diverging palette)
# Red for decrease, white for no change, green for increase
diff_vis = {
  'min': -0.2, # Adjust these based on expected magnitude of change
  'max': 0.2,  # A common range for NDVI change, can be adjusted for more contrast
  'palette': ['red', 'orange', 'yellow', 'white', 'lightgreen', 'green', 'darkgreen']
}


# --- Step 1: Compute Annual Median NDVI Images (Robust GEE Way) ---
# Note: In Python, you can often define a lambda or a nested function for map,
# or define a top-level function like apply_scale_factors. Here we use a regular function.
def compute_annual_median_ndvi(year_ee):
  year_ee = ee.Number(year_ee)

  year_string = year_ee.format('%d')

  end_date_current_year = ee.Algorithms.If(
    year_ee.eq(2025),
    year_string.cat('-06-30'),
    year_string.cat('-12-31')
  )

  start_date_current_year = year_string.cat('-01-01')

  year_collection = ndvi_collection \
    .filterDate(start_date_current_year, end_date_current_year)

  # DEBUG: Use server-side string concatenation for print statements within map
  # In Python, print inside a map function doesn't work directly like JS `print`.
  # You'd typically print sizes *after* mapping, or use .getInfo() if truly needed.
  # For demonstration, we'll keep the .getInfo() for direct debugging.
  # print(ee.String('Collection for ').cat(year_string).cat(' before median() size:').cat(year_collection.size().format()).getInfo()) # REMOVED .getInfo()

  median_ndvi = ee.Algorithms.If(
    year_collection.size().gt(0),
    year_collection.median().select('NDVI'),
    ee.Image(0).rename('NDVI').selfMask()
  )

  return ee.Image(median_ndvi).set('year', year_ee)

annual_median_ndvi_images_list = years_list.map(compute_annual_median_ndvi)


# NEW DEBUGGING: Evaluate the lists before creating the dictionary
# In Python, we use .getInfo() to evaluate. This is a synchronous call.
# Ensure keys are consistent strings (e.g., '2008', '2009', etc.)
keys = years_list.map(lambda year: ee.String(ee.Number(year).format('%d'))).getInfo()
if not keys:
    print('Error evaluating yearsList keys: Keys list is empty or evaluation failed.')
else:
    print('Client-side yearsList keys (strings):', keys) # Should contain '2025'

    images = annual_median_ndvi_images_list.getInfo()
    if not images:
        print('Error evaluating annualMedianNdviImagesList: Images list is empty or evaluation failed.')
    else:
        print('Client-side annualMedianNdviImagesList (contains image IDs/properties):', images)
        print('Length of client-side annualMedianNdviImagesList:', len(images))
        print('Expected length (endYear - startYear + 1):', end_year - start_year + 1)

        # Only proceed to create the dictionary and the next loop if evaluations succeed
        if len(keys) == len(images) and len(keys) == (end_year - start_year + 1):
            # Create ee.Dictionary from ee.List objects
            # Use the consistently formatted string keys
            annual_median_ndvi_images = ee.Dictionary.fromLists(
              years_list.map(lambda year: ee.String(ee.Number(year).format('%d'))), # Keys as ee.String (e.g., '2008')
              annual_median_ndvi_images_list # Values as ee.Image
            )

            print('Annual Median NDVI Images dictionary (server-side, GEE idiomatic):', annual_median_ndvi_images.getInfo())

            # NEW DEBUG: Check specific key after the dictionary creation, to ensure it's in the dictionary
            contains_end_year = annual_median_ndvi_images.contains(str(end_year)).getInfo()
            print(f'Is {end_year} key in dictionary?', contains_end_year)
            if contains_end_year:
                value_for_end_year = annual_median_ndvi_images.get(str(end_year)).getInfo()
                print(f'Value for {end_year} in dictionary:', value_for_end_year)

                # --- Step 2: Calculate Consecutive Year NDVI Differences (Server-Side) ---
                # This entire section is now nested inside the evaluation callbacks
                # This ensures `annual_median_ndvi_images` is fully resolved before used here.

                ndvi_difference_images_list = ee.List([])

                # Use a client-side for loop, as 'annual_median_ndvi_images' is now client-side
                for i in range(start_year, end_year): # Python range excludes the end_year
                    current_year1 = i
                    current_year2 = i + 1

                    # Cast to ee.Image as .get() returns an ee.Object
                    # Use consistent string keys
                    ndvi_year1 = ee.Image(annual_median_ndvi_images.get(str(current_year1)))
                    ndvi_year2 = ee.Image(annual_median_ndvi_images.get(str(current_year2)))

                    # Perform the subtraction.
                    ndvi_difference = ndvi_year2.subtract(ndvi_year1).rename('NDVI_Change')

                    # Add metadata (year1, year2) as properties to the image.
                    ndvi_difference = ndvi_difference.set({
                      'year1': current_year1,
                      'year2': current_year2,
                      'name': f'NDVI Change {current_year2}-{current_year1}'
                    })

                    # Add this annotated difference image to the list
                    ndvi_difference_images_list = ndvi_difference_images_list.add(ndvi_difference)

                # Evaluate the list of difference images once after the loop
                client_list_diff_images = ndvi_difference_images_list.getInfo() # Evaluate here
                if not client_list_diff_images:
                    print('*** List (ERROR) Details for Difference Images: ***')
                    print('Error evaluating NDVI difference images list: List is empty or evaluation failed.')
                    print('************************************************')
                else:
                    print('Client-side list of difference image properties (first 5 elements):', client_list_diff_images[:5])

                    # Create a Folium map centered on the AOI
                    map_center = aoi.centroid().coordinates().reverse().getInfo()
                    m = folium.Map(location=map_center, zoom_start=10)

                    # Add AOI to map
                    folium.GeoJson(aoi.getInfo(), name='Rimba Raya AOI', style_function=lambda x: {'color': 'red'}).add_to(m)


                    for image_properties in client_list_diff_images:
                        # Ensure image_properties and image_properties['id'] are valid before proceeding
                        # Access 'id' at the top level of the dictionary
                        if not image_properties or 'id' not in image_properties:
                            print('Skipping layer: Invalid image properties or ID for a difference image.')
                            continue

                        ee_image = ee.Image(image_properties['id'])
                        layer_name = image_properties['properties']['name'] # Access properties from dict

                        # Perform a reduceRegion to count non-masked pixels as a final check before adding layer
                        count_result = ee_image.reduceRegion(
                            reducer=ee.Reducer.count(),
                            geometry=aoi,
                            scale=30, # Use appropriate scale for your data (e.g., 30m for Landsat)
                            maxPixels=1e9 # Increase maxPixels if needed for large AOIs
                        ).get('NDVI_Change').getInfo()

                        if count_result and count_result > 0: # Check if count is not None and greater than 0
                            # Get the map ID and token from the Earth Engine image
                            map_id_dict = ee_image.getMapId(diff_vis)
                            tile_url = map_id_dict['tile_fetcher'].url_format

                            # Add the Earth Engine layer to the Folium map using TileLayer
                            folium.TileLayer(
                                tiles=tile_url,
                                attr='Map Data &copy; <a href="https://earthengine.google.com/">Google Earth Engine</a>',
                                name=layer_name,
                                overlay=True,
                                control=True
                            ).add_to(m)

                            print(f'Added layer: {layer_name} (Pixels: {count_result})')
                        else:
                            print(f'No data to display for: {layer_name} (No valid pixels or empty composite)')

                    # Add layer control to the map
                    folium.LayerControl().add_to(m)
                    # Display the map (if in a Jupyter/Colab environment)
                    display(m)

            else:
                print(f'Key {end_year} not found in dictionary.')
        else:
            print('Mismatch in list lengths. Dictionary will not be created correctly.')
            print('Years keys length:', len(keys))
            print('Images list length:', len(images))
            print('This typically means that `annualMedianNdviImagesList.map` failed for some years, resulting in fewer images than expected.')

# --- Optional: Add an overall Median NDVI map for the entire period as a baseline ---
overall_median_ndvi = ndvi_collection.median().select('NDVI').clip(aoi)

# Add overall median NDVI layer
# This converts the GEE image to a thumbnail URL for Folium
map_center = aoi.centroid().coordinates().reverse().getInfo()
m_overall = folium.Map(location=map_center, zoom_start=10)
folium.GeoJson(aoi.getInfo(), name='Rimba Raya AOI', style_function=lambda x: {'color': 'red'}).add_to(m_overall)

# Get the map ID and token for the overall median NDVI image
map_id_dict_overall = overall_median_ndvi.getMapId(ndvi_vis)
tile_url_overall = map_id_dict_overall['tile_fetcher'].url_format

# Add the overall median NDVI layer to the Folium map using TileLayer
folium.TileLayer(
    tiles=tile_url_overall,
    attr='Map Data &copy; <a href="https://earthengine.google.com/">Google Earth Engine</a>',
    name=f'Overall Median NDVI ({start_year}-{end_year})',
    overlay=True,
    control=True
).add_to(m_overall)

folium.LayerControl().add_to(m_overall)
display(m_overall)
tile_urls = {}
for y in range(start_year, end_year + 1):
    # Use consistent string key format
    y_key = str(y)
    # Get EE image for the year
    # Use consistent string key format
    ee_img = ee.Image(annual_median_ndvi_images.get(y_key))
    # Get tile URL for folium/ipyleaflet
    try:
        map_id = ee_img.getMapId(ndvi_vis)
        tile_urls[y] = map_id['tile_fetcher'].url_format
    except Exception as e:
        print(f"Could not get tile URL for year {y}: {e}")

# Step 2: Build the slider map
years = sorted(tile_urls.keys())
if years:
    center = map_center  # Already computed earlier
    m_slider = Map(center=center, zoom=10, scroll_wheel_zoom=True)

    # Prepare layers for each year
    tile_layers = {}
    for year in years:
        t = TileLayer(url=tile_urls[year], name=str(year), opacity=1 if year==years[0] else 0)
        tile_layers[year] = t
        m_slider.add_layer(t)

    m_slider.add_control(LayersControl(position='topright'))

    slider = widgets.IntSlider(
        value=years[0],
        min=years[0],
        max=years[-1],
        step=1,
        description='Year:',
        continuous_update=False
    )

    def on_year_change(change):
        for y in years:
            tile_layers[y].opacity = 1 if y == change['new'] else 0

    slider.observe(on_year_change, names='value')
    slider_control = WidgetControl(widget=slider, position='topright')
    m_slider.add_control(slider_control)

    display(m_slider)
else:
    print("No annual NDVI tile URLs available for slider map.")

print('Script complete. Check the map for NDVI change layers and Console for debug messages.')
