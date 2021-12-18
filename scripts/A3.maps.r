
## HOST SUITABILITY MAP
load(file="inputlyrs/rdata/land.rdata")
breaks <- c(0,20,40,60,80,100,999)
tags <- c("C10","C30", "C50", "C70", "C90", "OLD")
land$age.class <- cut(land$age, breaks=breaks, include.lowest=TRUE, right=TRUE, labels=tags)
land$sbw.suscep <- ifelse(land$spp %in% c("SAB", "OTH.RES.N") & land$age.class %in% c("C50", "C70", "C90", "OLD"), 1,
                          ifelse(land$spp %in% c("SAB", "OTH.RES.N") & land$age.class %in% c("C10", "C30"), 0.75,
                                 ifelse(land$spp %in% c("EPN", "OTH.RES.S") & land$age.class %in% c("C50", "C70", "C90", "OLD"), 0.5,
                                        ifelse(land$spp %in% c("EPN", "OTH.RES.S") & land$age.class %in% c("C10", "C30"), 0.25, 0))))
load(file="inputlyrs/rdata/mask.rdata")
map <- mask
map[!is.na(mask[])] <- land$sbw.suscep
plot(map, col=heat.colors(5)[5:1])




## VULNERABILITY MAP
source("R/sbw.vulnerability.r")
vuln <- vulnerability(land, MASK, 0.5)
MAP <- MASK
MAP[!is.na(MASK[])] <- vuln$v
# plot(MAP, col=c("grey80", heat.colors(5)[5:1]))
plot(MAP, col=heat.colors(5)[5:1])
group_by(vuln, spp, factor(host.pref)) %>% summarise(v=mean(v))