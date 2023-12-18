#---- Prep VPS Data for RSF Analysis ----
#elect only the columns needed by amt
dat <- grid_animalshort %>% dplyr::select(x='LON', y='LAT', t='DATETIME', id='TRANSMITTER')

#ensure all individuals have at least 100 positions
dat %>% 
  arrange(t) %>% 
  split(.$id) %>% 
  lapply(FUN=dim)

dat$t <- as.POSIXct(dat$t, format="%Y-%m-%d %H:%M:%S", tz="UTC")

dat_all <- dat %>% nest(-id)
dat_all <- dat_all %>% 
  mutate(trk=map(data, function(d) {
    amt::make_track(d, x, y, t, crs="EPSG:4326")
  }))

dat_all %>% 
  mutate(sr=lapply(trk, summarize_sampling_rate)) %>% 
  dplyr::select(id, sr) %>% 
  unnest(cols=c(sr))

dat1 <- dat_all %>% mutate(dat_clean=map(trk, ~ {
  .x %>% track_resample(rate=minutes(10), tolerance=seconds(120))
}))

#---- Reclassify substrate ----
#code provided in Hab.R

#---- RSF Analysis ----
dat_rsf <- dat1 %>% 
  mutate(rp=map(dat_clean, ~.x %>% random_points(n=nrow(.)*10) %>% 
                                   extract_covariates(substrate))) %>% 
  dplyr::select(id, rp) %>% unnest()

d <- rerun(100, dat_rsf %>% 
  dplyr::rename(sediment=5) %>% 
  mutate(fid=factor(id), sediment=factor(sediment)) %>% 
  mutate(sediment = fct_relevel(sediment, "1")) %>% #set same substrate reference level as used in iSSF
  mutate(i=1) %>% 
  group_by(id, case_) %>% 
  slice_sample(n=100) %>% 
  bam(case_~sediment+s(fid, bs="re"), 
      data=., 
      method="fREML", 
      family="binomial",
      discrete=T))

#create plot of relative selection strength
d <- d %>% 
  purrr::map(., ~broom::tidy(., parametric=T)) %>% 
  bind_rows(.) %>% 
  ggplot(aes(term, estimate %>% exp))+
  geom_hline(yintercept=1, lty=2)+
  geom_point(position="jitter", alpha=0.2)+
  geom_violin(fill="hotpink", colour=NA)+
  coord_flip()+
  scale_x_discrete(labels=c("Silt/mud with <50% gravel", "Coarse sediments", 
                            "Mixed sediments with veneer of mud", "Silt/mud","No data"))+
  labs(y="Relative Selection Strength", x="Substrate Class")+
  theme_classic()