### Script for extracting vegetation heights and lodging areas from 3D point clouds or digital elevation models ###
###################################################################################################################

# load packages
library(raster)
library(rgdal)
library(rgeos)
library(openxlsx)
library(lidR)



## get user input for processing settings ##

# create variable
clipping = " "
# path to the input data
wd = readline(prompt = "Enter the path to the input data (use / instead of \\): ")
# output path
output_path = readline(prompt = "Enter the output path (use / instead of \\): ")
# set wd
setwd(wd)
# name of the plot boundaries files
plot_boundaries_name = readline(prompt = "Enter the name of the plot boundary shapefile (without file extension): ")
# DTM or point cloud
input_type_DTM =  readline(prompt = "If you already have a DTM of the AOI enter d, if you have a point cloud enter p: ")
# name of the DTM or point cloud + desired resolution
if (input_type_DTM == "d"){
  DTM_name = readline(prompt = "Enter the name of the DTM (including the file extension): ")
}
if (input_type_DTM == "p"){
  DTM_points_name = readline(prompt = "Enter the name of the point cloud for the DTM (including the file extension): ")
  res_DTM = readline(prompt = "Enter the desired resolution for the DTM (in m): ")
  res_DTM = as.double(res_DTM)
  save_con_DTM = readline(prompt = "Save generated DTM (y/n)? ")
  if (save_con_DTM == "n"){
    save_con_DTM = FALSE
  } else{
    save_con_DTM = TRUE
  }
  clipping = readline(prompt = "Do you have a shapfile to clip the point cloud? If yes enter the name of the file (without file extension), if no enter n: ")
}
# DSM or point cloud
input_type_DSM =  readline(prompt = "If you already have a DSM of the AOI enter d, if you have a point cloud enter p: ")
# name of the DSM or point cloud + desired resolution
if (input_type_DSM == "d"){
  DSM_name = readline(prompt = "Enter the name of the DSM (including the file extension): ")
}
if (input_type_DSM == "p"){
  DSM_points_name = readline(prompt = "Enter the name of the point cloud for the DSM (including the file extension): ")
  res_DSM = readline(prompt = "Enter the desired resolution for the DSM (in m): ")
  res_DSM = as.double(res_DSM)
  save_con_DSM = readline(prompt = "Save generated DSM (y/n)? ")
  if (save_con_DSM == "n"){
    save_con_DSM = FALSE
  } else{
    save_con_DSM = TRUE
  }
  if (clipping == " "){
    clipping = readline(prompt = "Do you have a shapfile to clip the point cloud? If yes, enter the name of the file (without file extension), if no, enter n: ")
  }
}
# buffer size
buffer_setting = readline(prompt = "Set plot buffer width automatically (y/n)? ")
if (buffer_setting == "n"){
  buffer_width = readline(prompt = "Enter buffer width (in m): ")
}



## import data ##

# import/read DTM and or DSM if existing
if (input_type_DTM == "d"){
  DTM = raster(DTM_name)
  dot_pos = regexpr("\\.[^\\.]*$", DTM_name)
  DTM_name = substr(DTM_name, 1, dot_pos-1)
}

if (input_type_DSM == "d"){
  DSM = raster(DSM_name)
  dot_pos = regexpr("\\.[^\\.]*$", DSM_name)
  DSM_name = substr(DSM_name, 1, dot_pos-1)
}

# import plot boundaries
plot_boundaries = readOGR(".", plot_boundaries_name)




## create DEMs from 3D point clouds ##
 
# function to create DTM from point cloud
las_to_dtm = function(DTM_points, res_DTM, save_con_DTM, clip_file, DTM_name){
  if (mode(clip_file) != "character"){
    DEM <- clip_roi(DTM_points, clip_file)
    DEM <- classify_ground(DEM, algorithm = csf())
  } else{
    DEM <- classify_ground(DTM_points, algorithm = csf())
  }
  DEM <- grid_canopy(DEM, res = res_DTM, p2r())
  
  fill.na <- function(x, j=5) { if (is.na(x)[j]) { return(mean(x, na.rm = TRUE)) } else { return(x[j]) }}
  w <- matrix(1, 5, 5)
  
  DEM <- focal(DEM, w, fun = fill.na)
  
  if (save_con_DTM == TRUE){
    DTM_name = paste(DTM_name, "_dtm.tif", sep="")
    writeRaster(DEM, DTM_name, datatype ="FLT4S", overwrite=TRUE)
    print("DTM saved as raster")
  }
  else{print("DTM not saved as raster.")}
  
  return(DEM)
}

# function to create DSM from point cloud
las_to_dsm = function(DSM_points, res_DSM, DTM_z_min, DTM_Z_max, save_con_DSM, clip_file, DSM_name){
  DEM <- filter_poi(DSM_points, Z >= DTM_z_min, Z <= DTM_Z_max)
  if (mode(clip_file) != "character"){
    DEM <- clip_roi(DEM, clip_file)
  }
  DEM <- classify_ground(DEM, algorithm = csf())
  DEM <- grid_canopy(DEM, res = res_DSM, p2r())
  
  fill.na <- function(x, j=5) { if (is.na(x)[j]) { return(mean(x, na.rm = TRUE)) } else { return(x[j]) }}
  w <- matrix(1, 5, 5)
  
  DEM <- focal(DEM, w, fun = fill.na)
  
  if (save_con_DSM == TRUE){
    DSM_name = paste(DSM_name, "_dsm.tif", sep="")
    writeRaster(DEM, DSM_name, datatype ="FLT4S", overwrite=TRUE)
    print("DSM saved as raster")
  }
  else{print("DSM not saved as raster.")}
  
  return(DEM)
}

# import DTM and or DSM point cloud/s if existing 
if (input_type_DTM == "p"){
  DTM_points = readLAS(DTM_points_name)
  dot_pos = regexpr("\\.[^\\.]*$", DTM_points_name)
  DTM_name = substr(DTM_points_name, 1, dot_pos-1)
}
 
if (input_type_DSM == "p"){
  DSM_points = readLAS(DSM_points_name)
  dot_pos = regexpr("\\.[^\\.]*$", DSM_points_name)
  DSM_name = substr(DSM_points_name, 1, dot_pos-1)
  # get filter values 
  if (input_type_DTM == "D"){
    DTM_z_min = minValue(DTM) - 0.1
    DTM_z_max = maxValue(DTM) + 4
  } else{
    DTM_z_min = min(DTM_points@data[["Z"]]) - 0.1
    DTM_z_max = max(DTM_points@data[["Z"]]) + 4
  }
}

# import clip file if existing 
clip_file = " "
if (clipping != "n" & clipping != " "){
  clip_file = readOGR(".", clipping)
}

# set wd
wd = output_path
setwd(wd)

# create DTM and DSM raster from point cloud
if (input_type_DTM == "p" & input_type_DSM == "p"){
  DTM = las_to_dtm(DTM_points, res_DTM, save_con_DTM, clip_file, DTM_name)
  DSM = las_to_dsm(DSM_points, res_DSM, DTM_z_min, DTM_z_max, save_con_DTM, clip_file, DSM_name)
}

# disaggregate DTM if necessary
if (input_type_DTM == "p" & input_type_DSM == "p"){
  if (res_DSM < res_DTM){
    dis = res_DTM / res_DSM
    DTM = disaggregate(DTM, fact=dis, method='')
  }
}

# resample raster (with nearest neighbor) to align DTM with DSM, if necessary
if (isFALSE(compareRaster(DSM, DTM, stopiffalse = FALSE))) {
  DTM = projectRaster(DTM, DSM, method = "ngb")  
}

# Plot DTM and DSM
plot(DTM, main="DTM", legend.args=list(text="Elevation (m)", line=1.5), col=terrain.colors(100))
plot(DSM, main="DSM", legend.args=list(text="Elevation (m)", line=1.5), col=terrain.colors(100))



## calculate lodging ##

# create CHM
CHM = DSM-DTM

plot(CHM, main="CHM", legend.args=list(text="Height [m]", line=1.5), col=terrain.colors(100))


# clip CHM to plot-boundary extent and export it
clipped_CHM = mask(CHM, plot_boundaries)
clipped_CHM = crop(clipped_CHM, plot_boundaries)
CHM_name = paste(DSM_name, "_clipped__CHM.tif", sep="")
writeRaster(clipped_CHM, CHM_name, overwrite=TRUE)

# extract height values from CHM for each plot
extract_function = function(x, na.rm) c(percentile = quantile(x, .02, na.rm=TRUE), percentile = quantile(x, .98, na.rm=TRUE))
height_information = extract(clipped_CHM, plot_boundaries, fun = extract_function)

# convert height SPDF into df with values matching Plot-id
height_information = data.frame(id=plot_boundaries$id, height_information)

# merge info from Plots (e.g. Plot-id) with extracted height value table
height_information = as.data.frame(merge(plot_boundaries, height_information, by="id"))

# calculate difference between 98th and 2nd percentile 
height_information["diff_98_2"] = round(height_information$percentile.98. - height_information$percentile.2., digits = 2)

# add a column for the mode
height_information["modal_value"] = 0

# prepare variable and df needed in the loop
height_freq_df = data.frame(ID = 1:10000) 
height_freq_vec = c() 

# for loop to get mode
for (i in 1:nrow(height_information)){
  # clip plots
  c_c_CHM = mask(CHM, plot_boundaries[i, 1])
  c_c_CHM = crop(c_c_CHM, plot_boundaries[i, 1])
  # get frequency of the height values 
  frequency_m = freq(c_c_CHM, digits = 2)
  # convert matrix to data frame
  frequency_df = as.data.frame(frequency_m)
  # remove NA and NaN columns if necessary 
  remove_NA = 0
  if (is.nan(frequency_df[nrow(frequency_df),1]) == TRUE || is.nan(frequency_df[nrow(frequency_df)-1,1]) == TRUE){
    remove_NA = remove_NA + 1
  }
  if (is.na(frequency_df[nrow(frequency_df),1]) == TRUE || is.na(frequency_df[nrow(frequency_df)-1,1]) == TRUE){
    remove_NA = remove_NA + 1
  }
  frequency_df = frequency_df[1:(nrow(frequency_df)-remove_NA), 1:2]
  # get mode
  modal_value = frequency_df[frequency_df$count == max(frequency_df$count),1]
  if (length(modal_value) > 1){
    d = 0
    for (i2 in 1:(length(modal_value)-1)){
      d = modal_value[i2+1] - modal_value[i2] + d 
    }
    de = length(modal_value)*0.01-0.009
    if (d > de){
      modal_value = NA
     } else{
      modal_value = mean(modal_value)
    }
  }
  height_information[i, "modal_value"] = modal_value
}

# calculate difference between the percentiles and the modal value as well as the ratio between these differences
height_information["diff_mv_2"] = height_information$modal_value - height_information$percentile.2.
height_information["diff_98_mv"] = height_information$percentile.98. - height_information$modal_value
height_information["ratio"] = height_information$diff_mv_2 / height_information$diff_98_mv

# round a few columns 
height_information$percentile.2. = round(height_information$percentile.2., digits = 2)
height_information$percentile.98. = round(height_information$percentile.98., digits = 2)
height_information$diff_98_mv = round(height_information$diff_98_mv, digits = 2)
height_information$diff_mv_2 = round(height_information$diff_mv_2, digits = 2)
height_information$ratio = round(height_information$ratio, digits = 2)

# get max difference from the plots that have a ratio between the lowest ratio and the reciprocal value of the lowest ratio
max_diff = 0
min_ratio = min(height_information$ratio, na.rm = TRUE)
upper_ratio = min_ratio^-1
for (i in 1:nrow(height_information)){
  if (is.na(height_information[i, "ratio"]) == FALSE){
    if (height_information[i, "ratio"] >= min_ratio & height_information[i, "ratio"] <= upper_ratio & height_information[i, "diff_98_2"] > max_diff){
      max_diff = height_information[i, "diff_98_2"]  
    }
  }
}

# calculate threshold and upper and lower limit of uncertainty margin
height_information["threshold"] = round(height_information$percentile.98. - max_diff, digits = 2)
height_information["lower"] = round(height_information$threshold - abs(height_information$threshold)*0.05, digits = 2)
height_information["upper"] = round(height_information$threshold + abs(height_information$threshold)*0.05, digits = 2)

# calculation of the angle from which lodging is identified with high probability.
height_information["min_angle"] = round(acos(height_information$threshold / height_information$percentile.98.)*180/pi, digits = 2)

# for loop for the classification of each plot
for (i in 1:nrow(height_information)){
  # prepare reclass df
  c12 = height_information[i, "lower"]
  c23 = height_information[i, "threshold"]
  c34 = height_information[i, "upper"]
  c1 = 0
  if (c12<0){c1=c12}
  reclass_df = c(-1, 0, NA,
                 c1, c12, 1,
                 c12, c23, 2,
                 c23, c34, 3,
                 c34, Inf, 4)
  # clip plots
  c_c_CHM = mask(CHM, plot_boundaries[i, 1])
  c_c_CHM = crop(c_c_CHM, plot_boundaries)
  # classify plots and merge classified plots
  reclass_m = matrix(reclass_df,
                     ncol = 3,
                     byrow = TRUE)
  if (i == 1){
    CHM_classified = reclassify(c_c_CHM, reclass_m)
  } else{
    CHM_classified1 = reclassify(c_c_CHM, reclass_m)
    CHM_classified = merge(CHM_classified, CHM_classified1)
  }
}

# define extract function
extract_function = function(x, na.rm = TRUE) c(v_1 = sum(x == 1,  na.rm = TRUE), v_2 = sum(x == 2,  na.rm = TRUE), v_3 = sum(x == 3,  na.rm = TRUE), v_4 = sum(x == 4,  na.rm = TRUE))

# extract values
l_pp = extract(x = CHM_classified, y = plot_boundaries,df = TRUE, fun = extract_function)

# add number of raster cells column and plot id column  
l_pp["number of raster cells"] = c(1:nrow(l_pp))
l_pp["Plot"] = c(plot_boundaries$Parzelle)

# calculate number of raster cells
for (i in 1:nrow(l_pp)){
  l_pp[i, "number of raster cells"] = sum(l_pp[i, 4:ncol(l_pp)-2])
}

# reorder columns
l_pp = l_pp[, c(1, ncol(l_pp), ncol(l_pp)-1, 4:ncol(l_pp)-2)]

# change column names
colnames(l_pp)[4] = "most likely lodging"
colnames(l_pp)[5] = "probably lodging"
colnames(l_pp)[6] = "probably not lodging"
colnames(l_pp)[7] = "most likely not lodging"

# calculate class % per plot
for (i in 4:ncol(l_pp)){
  l_pp[paste(colnames(l_pp)[i], "%", sep=" ")] = round(l_pp[1:nrow(l_pp),i] / l_pp[1:nrow(l_pp),"number of raster cells"] * 100, 2)
}

# calculate class % for the whole CHM
class_df = data.frame(class = c("most likely lodging", "probably lodging", "probably not lodging", "most likely not lodging"))
n_raster_cells = sum(l_pp$`number of raster cells`)
for (i in 1:4){
  class_df[i, "%"] = round(sum(l_pp[1:nrow(l_pp), i+3]) / n_raster_cells * 100, digits = 2)
}

# create color lists for plotting
filling = c("red", "orange", "yellow", "green", "deeppink", "cyan", "blue", "darkorchid")
filling = filling[1:nrow(class_df)]

if(nrow(class_df) == 3) {
  color = c("red", "orange", "yellow", "green", "deeppink", "cyan", "blue", "darkorchid")
  color = color[1:nrow(class_df)]
} else {
  color = c("red", "orange", "yellow", "green", "deeppink", "cyan", "blue", "darkorchid")
  color = color[1:nrow(class_df)]
}

#add extra space to the right of the plot
par(mar=c(4, 3, 3, 7), xpd=TRUE)

# plot classified CHM
plot(CHM_classified,
     legend = FALSE,
     col = color,
     pch=1,
     main = "Classified CHM",
     xlab = "", 
     ylab = "")

legend("topright",
       inset=c(-0.4, 0),
       title="class",
       legend = class_df[1:nrow(class_df), c("class")],
       fill = filling,
       xpd=TRUE)

# export classifed_CHM
CHM_classified_name = paste(DSM_name, "_classified__CHM.tif", sep="")
writeRaster(CHM_classified, CHM_classified_name, overwrite=TRUE)



## extract height values from DSM ##

# buffer plot boundaries to exclude plot borders

# function for plot buffer with manual buffer width
manual_buffer = function(plot_boundaries, buffer_width){
  buffered_plots = vector("list", length(plot_boundaries))
  for (i in 1:length(plot_boundaries)) {
    plt = gBuffer(plot_boundaries[i,], width = buffer_width)
    cat("Plot", i, "buffered by", buffer_width, "m\n")
    plt$id = plot_boundaries[i,]$id
    plt$Parzelle = plot_boundaries[i,]$Parzelle
    buffered_plots[[i]] = plt
  }
  buffered_plots = do.call("rbind", buffered_plots)
  
  return(buffered_plots)
}

# function for plot buffer with automatic buffer width
# the loop tries different buffer widths (0.1, 0.3 and 0.5 m), starting from the highest
# if the resulting plot has an area more than 10% smaller than the original plot,
# the next lowest width is chosen. The minimum buffer is 0.1 m.
automatic_buffer = function(plot_boundaries){
  buffered_plots = vector("list", length(plot_boundaries))
  
  for (i in 1:length(plot_boundaries)) {
    buffer_test = gBuffer(plot_boundaries[i,], width = -0.5)
    
    if (area(buffer_test)/area(plot_boundaries[i,]) < 0.9){
      buffer_test = gBuffer(plot_boundaries[i,], width = -0.3)
      if (area(buffer_test)/area(plot_boundaries[i,]) < 0.9){
        plt = gBuffer(plot_boundaries[i,], width = -0.1)
        cat("Plot", i, "buffered by 10 cm\n")
        plt$id = plot_boundaries[i,]$id
        plt$Parzelle = plot_boundaries[i,]$Parzelle
        buffered_plots[[i]] = plt
      } 
      else {
        plt = gBuffer(plot_boundaries[i,], width = -0.3)
        cat("Plot", i, "buffered by 30 cm\n")
        plt$id = plot_boundaries[i,]$id
        plt$Parzelle = plot_boundaries[i,]$Parzelle
        buffered_plots[[i]] = plt}
    }
    else {
      plt = gBuffer(plot_boundaries[i,], width = -0.5)
      cat("Plot", i, "buffered by 50 cm\n")
      plt$id = plot_boundaries[i,]$id
      plt$Parzelle = plot_boundaries[i,]$Parzelle
      buffered_plots[[i]] = plt}
  }
  buffered_plots = do.call("rbind", buffered_plots)
  
  return(buffered_plots)
}

# buffer plots manually/automatically
if (buffer_setting == "y"){
  buffered_plots = automatic_buffer(plot_boundaries)
} else {
  buffer_width = as.numeric(buffer_width)
  buffered_plots = manual_buffer(plot_boundaries, buffer_width)
}

# overlay height and lodging rasters and se lodging areas to NA
masked_clipped_CHM = overlay(clipped_CHM, CHM_classified, fun = function(x,y){
  x[y<3] = NA
  return(x)
})

# extract height values from CHM for each plot
extract_function = function(x, na.rm) c(Mean_height = mean(x, na.rm=TRUE), 
                                        Median_height = median(x, na.rm=TRUE), 
                                        SD_height = sd(x, na.rm=TRUE), 
                                        MAD_height = mad(x, na.rm=TRUE), 
                                        percentile = quantile(x, .9, na.rm=TRUE))
extracted_heights = extract(masked_clipped_CHM, buffered_plots, fun = extract_function)

# convert height SPDF into df with values matching Plot-id
extracted_heights = data.frame(id=buffered_plots$id, extracted_heights)

# merge infos from Plots (e.g. Plot-id) with extracted height value table
extracted_heights = as.data.frame(merge(buffered_plots, extracted_heights, by="id"))

# convert height values to cm and print them
extracted_heights$Mean_height = extracted_heights$Mean_height * 100
extracted_heights$Median_height = extracted_heights$Median_height * 100
extracted_heights$SD_height = extracted_heights$SD_height * 100
extracted_heights$MAD_height = extracted_heights$MAD_height * 100
extracted_heights$percentile.90. = extracted_heights$percentile.90. * 100

# merge height values with lodging values in final table
colnames(extracted_heights)[2] = "Plot"
colnames(height_information)[2] = "Plot"
extracted_values_temp = merge(extracted_heights[,c("id", "Plot", "Mean_height", "Median_height", "SD_height", "MAD_height", "percentile.90.")],
                              l_pp[, c("Plot", "most likely lodging %", "probably lodging %", "probably not lodging %", "most likely not lodging %")],
                              by = "Plot", all.x=TRUE)
extracted_values = merge(extracted_values_temp, height_information[, c("Plot", "threshold", "lower", "upper", "min_angle")],
                         by = "Plot", all.x=TRUE)


# merge extract value df with plot SPDF and export it as SHP
extracted_values = extracted_values[,c(2,1,3:ncol(extracted_values))]
colnames(extracted_values)[1] = "ID"
colnames(extracted_values)[2] = "Parzelle"
buffered_plots = merge(plot_boundaries, 
                       extracted_values,
                       by = "Parzelle", all.x=TRUE)

writeOGR(plot_boundaries, "Plot Boundaries with extracted values", layer = "plot_boundaries_with_values",
         driver = "ESRI Shapefile", overwrite_layer = TRUE)

# export extracted values as xlsx
height_table_name = paste(DSM_name, "_Extracted_Values.xlsx", sep="")
write.xlsx(extracted_values, height_table_name, overwrite=TRUE)
print(extracted_values)

print("Finished processing without errors")
