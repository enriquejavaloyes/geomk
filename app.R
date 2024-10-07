#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinythemes)
library(shinycssloaders)
library(tidyverse)
library(sf)
library(osmdata)
library(mapSpain)
library(leaflet)
library(leaflet.extras)
library(spatstat)
library(pyramid)
library(sp)
library(spdep)
library(splancs)
library(dbscan)
library(tidygeocoder)
library(readxl)


sec_censal <- read_sf("data/SECC_CE_20240101.shp")

# Filtrar por municipio de Valencia
sec_valencia <- subset(sec_censal, CUMUN == "46250") %>%
  sf::st_transform(crs = 4326) %>%
  st_make_valid()

# Cargar datos demográficos y de renta
ambos_pob <- read_excel("data/edad_sexo_valencia.xlsx", sheet = "Ambos") %>%
  slice(-1) %>%
  rename(CUSEC = ...1) %>%
  mutate(CUSEC = substr(CUSEC, 1, 10)) %>%
  arrange(CUSEC)
pob_total_cens <- left_join(sec_valencia, ambos_pob, by = "CUSEC") %>%
  select(-c(2:18))

hombres_pob <- read_excel("data/edad_sexo_valencia.xlsx", sheet = "Hombres") %>%
  slice(-1) %>%
  rename(CUSEC = ...1) %>%
  mutate(CUSEC = substr(CUSEC, 1, 10)) %>%
  arrange(CUSEC)
pob_hombres_cens <- left_join(sec_valencia, hombres_pob, by = "CUSEC") %>%
  select(-c(2:18))

mujeres_pob <- read_excel("data/edad_sexo_valencia.xlsx", sheet = "Mujeres") %>%
  slice(-1) %>%
  rename(CUSEC = ...1) %>%
  mutate(CUSEC = substr(CUSEC, 1, 10)) %>%
  arrange(CUSEC)
pob_mujeres_cens <- left_join(sec_valencia, mujeres_pob, by = "CUSEC") %>%
  select(-c(2:18))

renta <- read_excel("data/ingreso_medio.xlsx") %>%
  rename(CUSEC = ...1) %>%
  mutate(CUSEC = substr(CUSEC, 1, 10)) %>%
  filter(grepl("^[0-9]+$", CUSEC)) %>%
  arrange(CUSEC)
renta_cens <- left_join(sec_valencia, renta, by = "CUSEC")

# Localizaciones
muni <- mapSpain::esp_get_munic() %>%
  sf::st_transform(crs = 4326)
valencia <- muni %>% filter(LAU_CODE == "46250")
bb_vlc <- sf::st_bbox(valencia)

# Datos de paradas de bus, metro y tranvía, gasolineras, y centros comerciales
bus_emt <- read_sf("data/emt.shp") %>%
  st_intersection(muni %>% filter(cmun == "250"))
coords <- sf::st_coordinates(bus_emt)
res <- dbscan(coords, eps = 0.0015, minPts = 1)
bus_emt$cluster <- data.frame(res$cluster)
bus_emt_clusters <- bus_emt %>%
  group_by(cluster) %>%
  summarise(geometry = st_centroid(st_combine(geometry)))

metro_bocas <- read_sf("data/fgv-bocas.shp") %>%
  st_intersection(muni %>% filter(cmun == "250"))
metro <- metro_bocas[grepl("[^46,]", metro_bocas$lineas), ]
coords <- sf::st_coordinates(metro)
res <- dbscan(coords, eps = 0.0015, minPts = 1)
metro$cluster <- data.frame(res$cluster)
metro_clusters <- metro %>%
  group_by(cluster) %>%
  summarise(geometry = st_centroid(st_combine(geometry)))

tranvia <- metro_bocas[grepl("4|6", metro_bocas$lineas), ]
coordst <- sf::st_coordinates(tranvia)
rest <- dbscan(coordst, eps = 0.0015, minPts = 1)
tranvia$cluster <- data.frame(rest$cluster)
tram_clusters <- tranvia %>%
  group_by(cluster) %>%
  summarise(geometry = st_centroid(st_combine(geometry)))


gasolineras_vlc <- bb_vlc %>%
  opq(timeout = 1000) %>%
  add_osm_feature(key = "amenity", value = "fuel") %>%
  osmdata_sf()
gasolineras_vlc_points <- gasolineras_vlc$osm_points %>%
  filter(!is.na(name)) %>%
  st_intersection(muni %>% filter(cmun == "250"))
gasolineras_vlc_points2 <- gasolineras_vlc$osm_polygons %>%
  st_intersection(muni %>% filter(cmun == "250")) %>%
  st_centroid()
gas_vlc <- rbind(
  gasolineras_vlc_points %>% select(name, addr.street, geometry),
  gasolineras_vlc_points2 %>% select(name, addr.street, geometry)
)

centros_comerciales_vlc <- bb_vlc %>%
  opq(timeout = 1000) %>%
  add_osm_feature(key = "shop", value = "mall") %>%
  osmdata_sf()
centros_comerciales_vlc_sf <- centros_comerciales_vlc$osm_polygons %>%
  st_intersection(muni %>% filter(cmun == "250"))
centros_comerciales_vlc_points <- st_centroid(centros_comerciales_vlc_sf)
corte_ingles <- read_excel("data/corte_ingles.xlsx")
corte_ingles_geocod <- tidygeocoder::geocode(corte_ingles, address = DIRECCION, method = "google") %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)

corte_ingles_geocod_sel <- corte_ingles_geocod %>% select(NOMBRE, DIRECCION, geometry)

centros_comerciales_vlc_points_sel <- centros_comerciales_vlc_points %>%
  rename(NOMBRE = name, DIRECCION = addr.street) 

centros_comerciales_vlc_points_sel <- centros_comerciales_vlc_points_sel %>% select(NOMBRE, DIRECCION, geometry)

cc_vlc <- rbind(corte_ingles_geocod_sel, centros_comerciales_vlc_points_sel) %>% 
  filter(NOMBRE != "El Manar")

cc_vlc_buffers <- st_buffer(cc_vlc, dist = 500)


# LISA --------------------------------------------------------------------

LISA <- function(data, variable, pesos, alpha=0.05) {
  P = localmoran(data[[variable]], listw = pesos)
  dif = data[[variable]] - mean(data[[variable]])
  lag = lag.listw(pesos, data[[variable]]) # Calcula el retardo (promedios)
  clag = dif - mean(lag) # Retardo - Media(Retardo)
  p = P[,5] # Toma la columna: Pr(z > 0) de P
  
  # Creamos un vector de n elementos (todos 5) y luego sustituimos
  quadrant = rep(5,length=nrow(P))+5
  quadrant[dif>0 & clag>0 & p<= alpha] = 1 # Alto-Alto
  quadrant[dif<0 & clag<0 & p<= alpha] = 2 # Bajo-Bajo
  quadrant[dif<0 & clag>0 & p<= alpha] = 3 # Bajo-Alto
  quadrant[dif>0 & clag<0 & p<= alpha] = 4 # Alto-Bajo
  
  # Grafico  
  brks = c(1,2,3,4,5)
  colors = c("red", "blue", "lightblue", "pink", "white")
  par(mar=rep(1,4))
  plot(data$geometry, border ="lightgray", col=colors[findInterval(quadrant, brks, all.inside=FALSE)], ylim = c(39.44, 39.505))
  legend("right", legend = c("High-High", "Low-Low", "Low-High", "High-Low", "Insignificant"),
         fill = colors, bty="n", cex=1.75, y.intersp=1, x.intersp=1)
}

# Aplicacion Shiny --------------------------------------------------------


ui <- fluidPage(
  theme = shinytheme("slate"),
  tags$head(
    tags$style(HTML("
      /* Estilo para la fuente League Spartan */
      @import url('https://fonts.googleapis.com/css?family=League+Spartan');
      body {
        font-family: 'League Spartan', sans-serif;
        font-size: 18px;
      }
    "))
  ),
  titlePanel("GEOMARKETING by Enrique Javaloyes"),
  sidebarLayout(
    sidebarPanel(
      selectInput("variable", "Seleccionar variable:", 
                  choices = c("Total Población" = "total_pob",
                              "Población Hombres" = "pob_hombres",
                              "Población Mujeres" = "pob_mujeres",
                              "Renta" = "renta")),
      width = 3,
      uiOutput("age_group_ui")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Mapa", 
                 leafletOutput("map", width = "100%", height = 700)),
        tabPanel("Información Adicional", 
                 fluidRow(
                   column(10, plotOutput("plot1")),
                   column(10, plotOutput("plot2"))
                 
                 ))
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Generar opciones de grupos de edad
  age_groups <- c("Total", "De 0 a 4 años", "De 5 a 9 años", "De 10 a 14 años", "De 15 a 19 años", 
                  "De 20 a 24 años", "De 25 a 29 años", "De 30 a 34 años", "De 35 a 39 años", 
                  "De 40 a 44 años", "De 45 a 49 años", "De 50 a 54 años", "De 55 a 59 años", 
                  "De 60 a 64 años", "De 65 a 69 años", "De 70 a 74 años", "De 75 a 79 años", 
                  "De 80 a 84 años", "De 85 a 89 años", "De 90 a 94 años", "De 95 a 99 años", 
                  "100 y más años")
  
  output$age_group_ui <- renderUI({
    if (input$variable != "renta") {
      selectInput("age_group", "Seleccionar grupo de edad:", choices = age_groups)
    }
  })
  
  filtered_data <- reactive({
    if (input$variable == "total_pob") {
      data <- pob_total_cens
    } else if (input$variable == "pob_hombres") {
      data <- pob_hombres_cens
    } else if (input$variable == "pob_mujeres") {
      data <- pob_mujeres_cens
    } else if (input$variable == "renta") {
      data <- renta_cens %>% select(CUSEC, geometry, last_col())
    }
    
    if (input$variable != "renta" && !is.null(input$age_group)) {
      age_group <- input$age_group
      if (age_group != "Total") {
        data <- data %>%
          select(CUSEC, geometry, all_of(age_group))
      }
    }
    data
  })
  
  output$map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles("CartoDB.Positron") %>%
      setView(lng = -0.37738226718744966, lat = 39.46954596430446, zoom = 13) %>%
      addCircleMarkers(data = bus_emt_clusters, color = "#FF0022", radius = 1, group = "Paradas de Bus") %>%
      addCircleMarkers(data = metro_clusters, color = "black", radius = 2, group = "Metro") %>%
      addCircleMarkers(data = tram_clusters, color = "#465E5E", radius = 2, group = "Tranvía") %>%
      addCircleMarkers(data = cc_vlc, color = "#00FFFF", group = "Centros Comerciales") %>% 
      addLayersControl(overlayGroups = c("Paradas de Bus" , "Metro", "Tranvía", "Gasolineras", "Centros Comerciales")) %>% 
      addMiniMap(tiles = providers$CartoDB.Positron,
                 position="bottomleft",
                 toggleDisplay = TRUE,
                 minimized = TRUE)
  })
  
  observe({
    req(filtered_data())
    data <- filtered_data()
    
    color_pal <- colorBin(c("#BE2A3E", "#F88F4D", "#90B960", "#22763F"), domain = as.numeric(data[[3]]), bins = 4) # COLOR SECCIONES
    
    leafletProxy("map") %>%
      clearShapes() %>%
      addCircles(data = gas_vlc, color = "#FF8200", radius = 30, group = "Gasolineras") %>%
      addPolygons(data = cc_vlc_buffers, color = "#00FFFF", fill = FALSE) %>%
      addPolygons(data = data, fillColor = ~color_pal(as.numeric(data[[3]])), fillOpacity = 0.25, weight = 1) %>%
      addPolygons(data = sec_valencia, fill = FALSE, color = "black", opacity = 0.3,weight = 1) %>%
      clearControls() %>% # Limpiar leyenda existente antes de añadir una nueva
      addLegend(
        position = "bottomright",
        colors = c("#FF0022", "black","#465E5E", "#FF8200", "#00FFFF"),
        labels = c("Paradas de Bus" , "Bocas de Metro", "Paradas de Tranvía", "Gasolineras", "Centros Comerciales"),
        title = "Leyenda",
        opacity = 0.5) %>%
      addLegend(
        position = "bottomright",
        pal = color_pal,
        values = as.numeric(data[[3]]),
        title = if (input$variable == "renta") "Renta Bruta media 2021" else "Población",
        opacity = 0.5)
  })
  
  output$plot1 <- renderPlot({
    req(filtered_data())
    data <- filtered_data()
    
    if (input$variable == "renta") {
      ggplot(data = renta_cens) +
        geom_sf(aes(fill = `Renta bruta media por persona 2021`)) +
        scale_fill_viridis_c(option = "plasma") +  # Escala de color, puedes elegir otras opciones como "viridis", "magma", etc.
        theme_void() +  # Tema minimalista para mejorar la apariencia
        labs( fill = "Renta Bruta") +
        coord_sf(ylim = c(39.44, 39.505), xlim = NULL)
      
    } else {
      variable <- input$age_group
      LISA(data, variable, pesos)
    }
                    
  })
  
  output$plot2 <- renderPlot({
    edades <- c("De 0 a 4 años", "De 5 a 9 años", "De 10 a 14 años", "De 15 a 19 años", 
                "De 20 a 24 años", "De 25 a 29 años", "De 30 a 34 años", "De 35 a 39 años", 
                "De 40 a 44 años", "De 45 a 49 años", "De 50 a 54 años", "De 55 a 59 años", 
                "De 60 a 64 años", "De 65 a 69 años", "De 70 a 74 años", "De 75 a 79 años", 
                "De 80 a 84 años", "De 85 a 89 años", "De 90 a 94 años", "De 95 a 99 años", 
                "100 y más años")
    
    hombres <- subset(hombres_pob, select = -c(1, 2))
    hombres <- apply(hombres, 2, sum)
    mujeres <- subset(mujeres_pob, select = -c(1, 2))
    mujeres <- apply(mujeres, 2, sum)
    
    
    datos<-data.frame(hombres,mujeres, edades)
    datos <- as.data.frame(datos[,-3])
    
    H1<-round(hombres/1000,0)
    M1<-round(mujeres/1000,0)
    datos<-data.frame(H1,M1,edades)
    pyramid(datos,Llab="Hombres",Rlab="Mujeres",Clab="Edad",main="Población Valencia 2021 \n (en miles)",Lcol="green", Rcol="cyan", Cgap=1, Cstep=2, Csize= 1.5)
  
  })
}

shinyApp(ui, server)





