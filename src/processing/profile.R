library(profvis)
library(bench)

profvis(readLog(files[1:2],useLocationList=TRUE,LocationList=LocationList,useFixes=FALSE))
profvis(makeSpatial(gnssData))
profvis(combineAsLines(gnssSF))
profvis(addAtmosphericDelay(gnssLines))
profvis(intersectBuilding(UCL,aboveGround))
