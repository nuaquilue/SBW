scn.name <- "TestSBW80"

spp.age <- read.table(paste0("outputs/", scn.name, "/SppByAgeClass.txt"), header=T)

spp.age.year <- mutate(spp.age, ymo=ifelse(age.class %in% c("C10", "C30", "C50"), "young",
                    ifelse(age.class %in% c("C70", "C90"), "mature", "old"))) %>% 
                group_by(run, year, spp, ymo) %>% summarise(area=sum(area)) %>% 
                group_by(year, spp, ymo) %>% summarise(area=mean(area)) 

sab.age.year <- filter(spp.age.year, spp=="SAB")

ggplot(spp.age.year, aes(x=year, y=area, group=ymo)) + geom_line(aes(colour=ymo)) +
  geom_point(aes(color=ymo)) + scale_color_brewer(palette="Dark2")  + theme_classic() +
  facet_wrap(~ spp)  

scn.name <- "test25"
track.sbw.defol.intens <- read.table(paste0("outputs/", scn.name, "/SBWdefoliation.txt"), header=T)

ggplot(filter(track.sbw.defol.intens, !is.na(pct)), aes(x=year, y=pct, group=as.factor(curr.intens.def))) +
  geom_line(aes(colour=as.factor(curr.intens.def))) + geom_point(aes(color=as.factor(curr.intens.def))) +
  scale_color_brewer(palette="Dark2")  + theme_classic()  +  facet_wrap(~ phase)  

ggplot(filter(track.sbw.defol.intens, !is.na(pct)), aes(x=year, y=ncell, group=as.factor(curr.intens.def))) +
  geom_line(aes(colour=as.factor(curr.intens.def))) + geom_point(aes(color=as.factor(curr.intens.def))) +
  scale_color_brewer(palette="Dark2")  + theme_classic()  #+  facet_wrap(~ phase)  
