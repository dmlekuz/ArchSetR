#rm(list=ls())

library(rcarbon)
library(data.table)
library(gstat)
library(automap)
library(rgdal)
library(ggplot2)
library(ggspatial)
library(GISTools)
library(scales)
#library(SDMTools)
library(zoo)
library(rgdal)
library(sp)
library(gridExtra)
library(raster)
library(spdep)
library(emstreeR)
library(dplyr) 

### variables
cnew<-"#c51b8a" #new grids color
cold<-"#fa9fb5" #old grids color
waterc<-"#9999AA" #water color
rollmeancolor<-"indianred"
fillcolor<-"gray80"

size <- 35000 #grid size
sizei<-10000 # interpolated grid size


###  prepare data
#set working directory
setwd("/Users/Dimitrij/Dropbox/neolitik/neo2")
## read dates and convert to dataframe
d <- read.table("alldates.csv", header=TRUE, sep=";")
dt<-data.table(d)
dt<-subset(dt, PERIOD   != "mesolithic")
dt<-dt[ ,c("ID","D","SD","SC")]
colnames(dt)<-c("rcode","d","sd","scode")
dt<-subset(dt, sd<=100) #only those with sd less then 100


## read shapefiles
sites<-readOGR(dsn="sites.shp",layer="sites") #sites
sea<-readOGR(dsn="sea.shp",layer="sea") #basemap
river<-readOGR(dsn="river.shp",layer="river") #basemap rivers
grid<-readOGR(dsn="grid.shp",layer="grid") #basemap grid
sa1<-readOGR(dsn="sa1.shp",layer="sa1") # large study area
#dem<-readGDAL(fname="dem1.tif")
#areas<-readOGR(dsn="areas.shp",layer="areas")
proj4string(sites)<-CRS("+init=epsg:4936")
proj4string(sea)<-CRS("+init=epsg:4936")
proj4string(river)<-CRS("+init=epsg:4936")
proj4string(grid)<-CRS("+init=epsg:4936")
proj4string(sa1)<-CRS("+init=epsg:4936")

#create grid
#hex_points <- spsample(sa1, type = "hexagonal", cellsize = size,n=1000, bb=10000000, offset=c(0.5,0.5))
#hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = size)
#hex_points <- spsample(sa1, type = "hexagonal", cellsize = size,n=1000, offset=c(0.2,0.2))
#hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = size)
#hex_ids<-as.data.frame(getSpPPolygonsIDSlots(hex_grid))
#grid<-SpatialPolygonsDataFrame(hex_grid, data=hex_ids,match.ID=F)
#colnames(grid@data)<-"id" # column ID contains grid reference
#hex_points$id<-over(hex_points, grid)


### Figure 1, 
### study area shaded dem with sites, rivers  and grid
###
waterc1<-"#CCCCDD" #uses different water colour for easthetic purposes

#ddd<-as.data.frame(dem,xy=T)

ggplot()+
  #geom_raster(data=ddd,aes(x=x,y=y, fill=band1),alpha=0.6)+
  #scale_fill_gradient(low="black",high="white",limits=c(20,1250),name=NULL,oob=squish)+
  layer_spatial(data=grid, alpha=0.2,fill="#aaaaaa",col=NA)+
  layer_spatial(data=sea,fill=waterc1,col=NA) +
  layer_spatial(data=river,col=waterc1, alpha=0.7) +
  layer_spatial(data=sites,size=0.3,col="black")+
  coord_sf(xlim=c(4320000,5250000),ylim=c(1890000,2650000))+
  theme(legend.position="none")+
  labs()+xlab("")+ylab("")

#ggsave("fig1.pdf",plot=last_plot())


### Figure 2
### Number of dates per grid
###  
# join dates with sites, grid
s<-data.frame(sites@coords,sites@data$sc) #convert sites to data frame
ds<-merge(dt,s, by.x="scode", by.y="sites.data.sc") #actual join
colnames(ds)<-c("scode", "rcode","d","sd","x","y")
#convert to spatialdataframe
coordinates(ds)<-cbind(ds$x,ds$y)
proj4string(ds)=CRS(proj4string(sites))
#add information in which grid is date
ds@data$grid<-over(ds, grid[, 'id'])$id #get ids of grids date is in
ds@data$region<-over(ds, grid[, 'id'])$region

#only sites that have dates
sites<-subset(sites, sites$sc %in% unique(ds$scode))


gr.dat<-aggregate(ds@data$grid,by=list(ds@data$grid),FUN=length)
colnames(gr.dat)<-c("id","count")
grid.d<-merge(grid,gr.dat, by.x="id", by.y="id")
grid.d<-subset(grid.d, grid.d@data$count>0)

ggplot()+
  layer_spatial(data=grid, alpha=0.2,fill="#aaaaaa",col=NA)+
  layer_spatial(data=grid.d,col=NA,aes(fill=count))+
  scale_fill_distiller(palette="RdPu", direction=1,oob=squish,trans="log10") +
  layer_spatial(data=sea,fill=waterc,col=NA,alpha=0.8) +
  layer_spatial(data=river,col=waterc, alpha=0.8) +
  layer_spatial(data=sites,size=0.2)+
  coord_sf(xlim=c(4320000,5250000),ylim=c(1890000,2650000))+
  labs(title="",fill="N")+xlab("")+ylab("")

#ggsave("fig2.pdf",plot=last_plot()) 



### Fig 3 ALT
### earliest apperance in a grid using SPD
###


binSize<-100 #bin size used in spds
cutoff<-0.05 #2 sigma

df<-as.data.table(ds)
gridss<-df[,.(.N),by=grid]$grid #list of grids
#this will be our output
out<-data.frame(grid=character(),d=integer(),num=integer(),stringsAsFactors = F ) #num is number of dates

# calculate spds for each grid
for(g in gridss) {
          print(g)
     ss<-subset(df,df$grid==g) #select all dates from a particular grid g
     ss.cd=calibrate(x=ss$d,errors=ss$sd,calCurves='intcal20',ncores=3) #calibrate them
     ss.bins = binPrep(sites=ss$grid,ages=ss$d,h=binSize, method = "complete") # create bins
     ss.spd.bins = spd(timeRang=c(10000,3000),ss.cd,bins=ss.bins,spdnormalised=T)
     #now find edge of lower 5% (lower 2sigma) 
     #accumulate probabilities until it reaches 2.5%
     s<-0
     print(g)
     for(n in seq(from=1,to=nrow(ss.spd.bins$grid))) {
       s<-s+ss.spd.bins$grid[n,2]
       if(s>cutoff/2) { 
         dte<-ss.spd.bins$grid[n ,1]
         break }
     }
     # aa<-subset(..ss.spd.bins$grid, ss.spd.bins$grid$PrDens > probTwoSigma)
     out<-rbind(out,data.frame(grid=g, d=dte, num=nrow(ss)))

 } 
 
#get spatial extent of grids
grid.d<-merge(grid,out, by.x="id", by.y="grid")
grid.d<-subset(grid.d, grid.d@data$d>0) #select only grid cells with data

ggplot()+
  layer_spatial(data=grid, alpha=0.2,fill="#aaaaaa",col=NA)+
  layer_spatial(data=grid.d,col=NA,aes(fill=d))+
  scale_fill_distiller(palette="YlOrBr",direction=1,limits=c(6500,8100),oob=squish) +
  layer_spatial(data=sea,fill=waterc,col=NA,alpha=0.8) +
  layer_spatial(data=river,col=waterc, alpha=0.8) +
  layer_spatial(data=sites,size=0.2)+
  coord_sf(xlim=c(4320000,5250000),ylim=c(1890000,2650000))+
  labs(title="",fill="cal BP")+xlab("")+ylab("")+ annotation_scale()

#ggsave("mlekuz_slika3.pdf",plot=last_plot()) 


### Fig 4 
### krigging
###
points <- SpatialPointsDataFrame(coords=grid, data=grid@data, proj4string=CRS("+init=epsg:4936"))
colnames(points@data)<-"id"
points.d<-merge(points,grid.d, by="id") #centroids of grids with the earliest date
pts.all<-data.frame(points.d)
pts<-subset(pts.all,d>0 & num>1) #dataframe with centroids of points  that have earliest date AND have more than one date
# convert back to spatialdataframe
coordinates(pts)<-~ coords.x1 + coords.x2

#create grid of interpolated points over study area
p.i <- spsample(sa1, type = "regular", cellsize = sizei,n=1000) #type = "hexagonal"
#coordinates(p.i)<- ~ x + y
proj4string(p.i)=CRS(proj4string(pts))

#data.kriged<-autoKrige(d~1,pts,p.i,kappa=0.88,model="Mat",fix.values=c(0,125000,70000))
#model = c("Sph", "Exp", "Gau", "Ste")
data.kriged<-autoKrige(d~1, pts, p.i, fix.values=c(300000,90000,100000), model="Exp")
plot(data.kriged)
kriged<-as.data.frame(data.kriged$krige_output)
names(kriged)<-c("x1","x2","d.pred","d.var","var.stdev")
#kriged<-subset(kriged,var.stdev<850)

ggplot()+
  layer_spatial(data=grid, alpha=0.2,fill="#aaaaaa",col=NA)+
  geom_contour(data=kriged,aes(x1,x2,z=d.pred,colour=stat(level)),binwidth=100, size=0.66)  +
  scale_colour_distiller(palette="YlOrBr",direction=1, limits=c(7000,8000),oob=squish) +
  layer_spatial(data=sea,fill=waterc,col=NA,alpha=0.7) +
  layer_spatial(data=river,col=waterc, alpha=0.8) +
  layer_spatial(data=sites,size=0.2)+
  coord_sf(xlim=c(4320000,5250000),ylim=c(1890000,2650000))+
  labs(title="",colour="cal BP")+xlab("")+ylab("")+
        annotation_scale()

#ggsave("mlekuz_slika4.pdf",plot=last_plot())   




### SPDs
### 
###     

grids.spd<-NULL

for(g in gridss){
     ss<-subset(df,df$grid==g) #select all dates from a particular grid g
     ss.cd=calibrate(x=ss$d,errors=ss$sd,calCurves='intcal20',ncores=3) #calibrate them
     ss.bins = binPrep(sites=ss$grid,ages=ss$d,h=binSize) # create bins
     ss.spd.bins = spd(ss.cd,timeRange=c(11000,2500),bins=ss.bins,spdnormalised=T)
     #get 95% cutoff value
     #now find edge of lower 5% (lower 2sigma) 
     #accumulate probabilities until it reaches 2.5%
     s<-0
     val<-0 # density value which determines 2sigma cutoff #we need to find it
     for(n in seq(1,nrow(ss.spd.bins$grid))) {
       s<-s+ss.spd.bins$grid[n,2]
       if(s>cutoff/2) { 
         val<-ss.spd.bins$grid[n,2]
         break 
       }
     }   
     print(g)
     print(val)
     ss.spd.bins$id<-g
     ss.spd.bins$cutoff<-val
     grids.spd<-rbind(grids.spd,ss.spd.bins)
}



grids.at<-function(date,agrids)
{
  gridsat<-NULL
  for(n in 1:dim(agrids)[1])
  {
   g=agrids[n,]  
   if(g$grid$PrDens[g$grid$calBP==date]>g$cutoff)
     { gridsat<-c(gridsat,g$id)}
  }  
  return(gridsat)  
}

g<-grids.at(5600,grids.spd)
g1<-c(g1,g)
gds<-subset(grid, grid$id %in% g1)

ggplot()+
  layer_spatial(data=grid, alpha=0.2,fill="#aaaaaa",col=NA)+
  layer_spatial(data=gds,fill="#aa2200",col=NA)+
  scale_fill_distiller(palette="YlGnBu",direction=1,limits=c(6500,8500),oob=squish) +
  layer_spatial(data=sea,fill=waterc,col=NA,alpha=0.8) +
  layer_spatial(data=river,col=waterc, alpha=0.8) +
  layer_spatial(data=sites,size=0.2)+
  coord_sf(xlim=c(4320000,5250000),ylim=c(1890000,2650000))+
  labs(title="",fill="cal BP")+xlab("")+ylab("")+ annotation_scale()


  ggsave("5600.pdf",plot=last_plot())     


### Fig 5 spd plot
### 
###  

get_density<-function(reg)
# needs grid as  input
{
  gridss<-subset(grid, region==reg)$id #list of grids
  dens.all<-0

  for(g in gridss) {
     ss<-subset(df,df$grid==g) #select all dates from a particular grid g
     ss.cd=calibrate(x=ss$d,errors=ss$sd,calCurves='intcal20',ncores=3) #calibrate them
     ss.bins = binPrep(sites=ss$grid,ages=ss$d,h=binSize) # create bins
     ss.spd.bins = spd(ss.cd,timeRange=c(10000,3500),bins=ss.bins,spdnormalised=T)
     print(g)
     dens.all<-dens.all+ss.spd.bins$grid$PrDens
  }
 #dens<-data.frame(ss.spd.bins$grid$calBP,dens.all)
 #colnames(dens)<-c("calBP","all")
 return(dens)

}





dens.all<-0
for(g in gridss) {
     ss<-subset(df,df$grid==g) #select all dates from a particular grid g
     ss.cd=calibrate(x=ss$d,errors=ss$sd,calCurves='intcal20',ncores=3) #calibrate them
     ss.bins = binPrep(sites=ss$grid,ages=ss$d,h=binSize) # create bins
     ss.spd.bins = spd(ss.cd,timeRange=c(10000,3500),bins=ss.bins,spdnormalised=T)
     print(g)
     dens.all<-dens.all+ss.spd.bins$grid$PrDens
}



dens<-data.frame(ss.spd.bins$grid$calBP,dens.all)
colnames(dens)<-c("calBP","all")

ggplot()+
   geom_area(data=dens,aes(x=calBP,y=all))+
   scale_x_reverse()+
     geom_line(data=dens, aes(x=calBP,y=rollmean(all,200,na.pad=TRUE),colour="#aa2200"), linetype=1,size=1)+
   scale_x_reverse()+
   #geom_line(data=dens, aes(x=dens$calBP,y=rollmean(dens$all,200,na.pad=TRUE),colour="grey"), linetype=1,size=1)+
   xlim(8500,5000)+
   labs(title="")+xlab("cal BP")+ylab("Summed probability")

ggsave("summedProbNeol.pdf",plot=last_plot())     


### Minimum spanning network
### 
###  

get.point.id<-function(x,y, points)
{
  for(i in 1:nrow(points)){
    if(points[i,]$x==x | points[i,]$y==y){
      break()
    } 
   }
    return(points[i,]$id)   
}


points <- SpatialPointsDataFrame(coords=grid, data=grid@data, proj4string=CRS("+init=epsg:4936"))
#points<-hex_points
grf.all<-data.frame("", "", 0, 0, 0) #all minimum soaning trees together
colnames(grf.all)<-c("f","t","d","c","y") #from,to,distance,count, year 

for(year in seq(8500,5000,-1))
{
  g.at<-data.frame(grids.at(year,grids.spd)) #grids occupied at year
  if (nrow(g.at)>1)
  {
    colnames(g.at)<-"id"
    gat<-subset(points,points$id %in% g.at$id)# centorid of grids occupied at year


    coords<-coordinates(gat)
    #black magic -- we add some outside celles, because ComputeMST does not
    # work with les then few nodes
      for(i in 1:10)
      {
          coords<-rbind(coords, c(5141545+100000*i,1924742+100000*i))
      }

    mst<-ComputeMST(coords) #minimum spanning tree at year

    #we have to repair mst because ComputeMst sometimes reverses coordinates
    #for(j in 1:nrow(mst))
    #{ 
    #  if(mst[j,]$gat.x<mst[j,]$gat.y)
    #  { 
    #    a<-mst[j,]$gat.x
    #    mst[j,]$gat.x<-mst[j,]$gat.y
    #    mst[j,]$gat.y<-a
    #  }
    #}  
    
    
    mst.p<-SpatialPointsDataFrame(coords=as.matrix(mst[,c(1,2)]),data=data.frame(mst[,c(3,4,5)])) #identify grids of minimum spaning tree nodes 
    projection(mst.p)<-CRS(proj4string(points))
    mst.p$grd<-over(mst.p,grid)$id
   

    m<-data.frame(mst.p)

    grf<-data.frame(m[m$from,]$grd, m[m$to,]$grd,m$distance) #convert references to grid ids
    colnames(grf)<-c("t","f","d")
    grf <- grf[-nrow(grf),] # delete last one whic is irrelevant
    grf$c<-0
    grf$y<-year #add current year


    grf.all<-rbind(grf.all, grf) #add to all minumum spanning treees
  }
}  


#filter galore!!!!
grf.all <- grf.all[-1,] #remove first because it is not ok
#grf.all <- grf.all %>% distinct(f,t,.keep_all=T) #only distinct
#grf.all <- subset(grf.all, grf.all$f!=NA & grf.all$t!=NA)
grf.all<-grf.all %>% filter(grf.all$t!="NA" & grf.all$f!="NA" ) #remove NAs


a<-grf.all %>% group_by(t,f) %>% filter(y == max(y)) #find oldest

b<-grf.all<-grf.all %>% count(t, f, d) #count number of connection

#merge count data and oldest data, update col names
grf.all<-merge(a, b, by.x=c("t","f"), by.y=c("t","f"))
grf.all<-grf.all[c("t", "f", "d.x","y","n")]
colnames(grf.all)<-c("t", "f", "d","y","n")

#only those with large n of connections
#grf.all<-subset(grf.all,n>10)
#only short connections, les than 500 km

grf.all<-subset(grf.all,d<400000)

# only older than 7500
#grf.all<-subset(grf.all,y>8000)

l1<-SpatialLinesDataFrame(sl=SpatialLines(LinesList=list(Lines(list(Line(coords=matrix(c(0,1,0,1),2,2))),ID=1))),
    data=data.frame(f="",t="",d=0,c=0,y=0))[-1,] #empty
proj4string(l1)<-CRS("+init=epsg:4936")



for(i in 1:nrow(grf.all))
{  
  print(grf.all[i,]$f)
  f<-points[points$id %in% grf.all[i,]$f,]
  t<-points[points$id %in% grf.all[i,]$t,]
  l<-as(rbind(f,t),"SpatialLines")

  l1<-rbind(l1,SpatialLinesDataFrame(l,data=as.data.frame(grf.all[i,]),match.ID=F))

}  



ggplot()+
  layer_spatial(data=l1,aes(size=n,color=y)) +
  scale_size(range = c(0.0, 1.5))+
  scale_colour_distiller(palette="YlOrBr",direction=1, limits=c(6500,8500),oob=squish) +
  layer_spatial(data=grid, alpha=0.2,fill="#aaaaaa",col=NA)+
  layer_spatial(data=sea,fill=waterc,col=NA,alpha=0.7) +
  layer_spatial(data=river,col=waterc, alpha=0.8) +
  layer_spatial(data=sites,size=0.2)+
  coord_sf(xlim=c(4320000,5250000),ylim=c(1890000,2650000))+
  labs(title="")+xlab("")+ylab("")

     ggsave("MST.pdf",plot=last_plot())   

### mInimal spanning tree analysis
###
###


grf.all<-data.frame("", "", 0, 0, 0) #all minimum soaning trees together
colnames(grf.all)<-c("f","t","d","c","y") #from,to,distance,count, year 

for(year in seq(8500,5000,-1))
{
  g.at<-data.frame(grids.at(year,grids.spd)) #grids occupied at year
  if (nrow(g.at)>1)
  {
    colnames(g.at)<-"id"
    gat<-subset(points,points$id %in% g.at$id)# centorid of grids occupied at year


    coords<-coordinates(gat)
    #black magic -- we add some outside celles, because ComputeMST does not
    # work with les then few nodes
      for(i in 1:10)
      {
          coords<-rbind(coords, c(5141545+100000*i,1924742+100000*i))
      }

    mst<-ComputeMST(coords) #minimum spanning tree at year

    #we have to repair mst because ComputeMst sometimes reverses coordinates
    #for(j in 1:nrow(mst))
    #{ 
    #  if(mst[j,]$gat.x<mst[j,]$gat.y)
    #  { 
    #    a<-mst[j,]$gat.x
    #    mst[j,]$gat.x<-mst[j,]$gat.y
    #    mst[j,]$gat.y<-a
    #  }
    #}  
    
    
    mst.p<-SpatialPointsDataFrame(coords=as.matrix(mst[,c(1,2)]),data=data.frame(mst[,c(3,4,5)])) #identify grids of minimum spaning tree nodes 
    projection(mst.p)<-CRS(proj4string(points))
    mst.p$grd<-over(mst.p,grid)$id
   

    m<-data.frame(mst.p)

    grf<-data.frame(m[m$from,]$grd, m[m$to,]$grd,m$distance) #convert references to grid ids
    colnames(grf)<-c("t","f","d")
    grf <- grf[-nrow(grf),] # delete last one whic is irrelevant
    grf$c<-0
    grf$y<-year #add current year


    grf.all<-rbind(grf.all, grf) #add to all minumum spanning treees
  }
}

#filter galore!!!!
grf.all <- grf.all[-1,] #remove first because it is not ok
grf.all <- grf.all %>% distinct(f,t,.keep_all=T) #only distinct
#grf.all <- subset(grf.all, grf.all$f!=NA & grf.all$t!=NA)
grf.all<-grf.all %>% filter(grf.all$t!="NA" & grf.all$f!="NA" ) #remove NAs
#grf.all<-subset(grf.all,d<400000)



grf.all <- grf.all %>% distinct(f,.keep_all=T) #this is key, we keep only expansions from allready established
 
# lenghth of spanning network by  year
grf.l<-grf.all %>% group_by(y) %>% summarise(length=max(d),n=n())
grf.l$length<-grf.l$length/1000 #in km

grf.l<-arrange(grf.l,desc(y))


ggplot()+
   geom_line(data=grf.l,aes(x=y,y=cumsum(length)))+
   #geom_line(data=grf.all, aes(x=y,y=rollmean(d,5,na.pad=TRUE),colour="#aa2200"), linetype=1,size=1)+
   scale_x_reverse()+
   #geom_line(data=dens, aes(x=dens$calBP,y=rollmean(dens$all,200,na.pad=TRUE),colour="grey"), linetype=1,size=1)+
   xlim(8500,5000)+
   labs(title="")+xlab("cal BP")+ylab("Expansion length")

   ggsave("MSTLength1.pdf",plot=last_plot())    








grf.m<-grf.all %>% group_by(y) %>% summarise(maxlen=max(d),n=n())

ggplot()+
   geom_area(data=grf.m,aes(x=y,y=maxlen))+
   geom_line(data=grf.m, aes(x=y,y=rollmean(maxlen,20,na.pad=TRUE),colour="#aa2200"), linetype=1,size=1)+
   scale_x_reverse()+
   #geom_line(data=dens, aes(x=dens$calBP,y=rollmean(dens$all,200,na.pad=TRUE),colour="grey"), linetype=1,size=1)+
   xlim(8500,5000)+
   labs(title="")+xlab("cal BP")+ylab("MST length")

      ggsave("MSTexpansion.pdf",plot=last_plot())  



