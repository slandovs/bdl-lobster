#---- Plot tracks on substrate raster ----
labs <- expression("Coarse sediments", "Mixed sediments with veneer of mud", "Silt/mud with <50% gravel",
                   "Silt/mud", "No data")
hab_dat <- p +
  coord_fixed(xlim=c(min(grid_animalshort$LON), max(grid_animalshort$LON)),
              ylim=c(min(grid_animalshort$LAT), max(grid_animalshort$LAT))) +
  geom_point(data = grid_animalshort,
             aes(LON, LAT, colour=TRANSMITTER), size=1.5, inherit.aes = F) +
  scale_fill_manual(labels=labs, values=c("#D6BEA1","#7AB6F5","#AACCE3","#FFFFCE", "gray90")) + 
  labs(x="Longitude", y="Latitude", colour="Lobster ID") +
  theme(plot.title = element_text(hjust = 0.5, vjust=3), legend.text.align = 0) 

#---- Choose a sampling rate and tolerance ----
#first, resample all data to form regular bursts: this is required for SSFs because selection is not scale invariant (Signer et al 2019)
#thus, sampling rates should be similar for different animals in a given study (use the highest median as your rate)

dat <- grid_animalshort %>% dplyr::select(x='LON', y='LAT', t='DATETIME', id='TRANSMITTER')
dat_all <- dat %>% nest(-id)

dat_all <- dat_all %>% 
  mutate(trk=map(data, function(d) {
    amt::make_track(d, x, y, t, crs=sp::CRS("+init=epsg:4326"))
  }))
dat_all %>% mutate(sr=lapply(trk, summarize_sampling_rate)) %>% dplyr::select(id, sr) %>% unnest(cols=c(sr))
#highest median is 4, highest q3 is 10

#---- Reclassify substrate ----
substrate[is.na(substrate[])] <- 5

reclass_substrate <- function(s) {
  fct_collapse(factor(s),
               till = "0",
               mudL = "1",
               csed = "2",
               ssed = c("3","4"),
               nd = "5")
}

#---- Functions for iSSFs and bootstrapping ----
issf <- function(i){
  t <- amt::make_track(i, .x=LON, .y=LAT, .t=DATETIME, crs=sp::CRS("+init=epsg:4326")) 
  r <- track_resample(t, rate=minutes(10), tolerance=minutes(2)) %>% 
    filter_min_n_burst(min_n=3) %>% 
    steps_by_burst() %>% 
    random_steps(n_control = 100) %>%
    extract_covariates(substrate, where="both") %>% 
    mutate(substrate_start = reclass_substrate(layer_start),
           substrate_end = reclass_substrate(layer_end),
           cos_ta_ = cos(ta_),
           log_sl_ = log(sl_)) %>%
    filter(!is.infinite(log_sl_)) %>%
    filter(!is.na(ta_))
}

boot <- function(i){
  b <- i %>% bootstrap(., 1000)
  a <- b %>% .$strap %>% 
    purrr::map( ~ fit_issf(case_ ~ substrate_end + log_sl_ + cos_ta_ + 
                             substrate_start:log_sl_ + substrate_start:cos_ta_ +
                             strata(step_id_), data = .)) %>% 
    map(., first) %>% 
    purrr::map(~broom::tidy(.)) %>% 
    bind_rows()
}

#---- issf function ----
issf27 <- function(i){
  t <- amt::make_track(i, .x=LON, .y=LAT, .t=DATETIME, crs=sp::CRS("+init=epsg:4326")) 
  r <- track_resample(t, rate=minutes(10), tolerance=minutes(2)) %>% 
    filter_min_n_burst(min_n=3) %>% 
    steps_by_burst() %>% 
    time_of_day(include.crepuscule=FALSE) %>% 
    random_steps(n_control = 100) %>% 
    extract_covariates(substrate, where="both") %>%
    filter(layer_end != "5") %>% 
    filter(layer_start != "5") %>% 
    mutate(substrate_start = reclass_substrate(layer_start),
           substrate_end = reclass_substrate(layer_end),
           cos_ta_ = cos(ta_), 
           log_sl_ = log(sl_)) %>%
    filter(!is.infinite(log_sl_)) %>% 
    filter(!is.na(ta_))
}
issf027 <- issf27(grid_animalshort %>% filter(TRANSMITTER=="BL 027"))

issf030 <- issf(grid_animalshort %>% filter(TRANSMITTER=="BL 030"))

issf32 <- function(i){
  t <- amt::make_track(i, .x=LON, .y=LAT, .t=DATETIME, crs=sp::CRS("+init=epsg:4326")) 
  r <- track_resample(t, rate=minutes(10), tolerance=minutes(2)) %>% 
    filter_min_n_burst(min_n=3) %>% 
    steps_by_burst() %>% 
    time_of_day(include.crepuscule=FALSE) %>% 
    random_steps(n_control = 100) %>% 
    extract_covariates(substrate, where="both") %>% 
    filter(layer_end != "2") %>% 
    filter(layer_start != "2") %>% 
    mutate(substrate_start = reclass_substrate(layer_start),
           substrate_end = reclass_substrate(layer_end),
           cos_ta_ = cos(ta_), 
           log_sl_ = log(sl_)) %>% 
    filter(!is.infinite(log_sl_)) %>% 
    filter(!is.na(ta_))
}
issf032 <- issf(grid_animalshort %>% filter(TRANSMITTER=="BL 032"))

issf34 <- function(i){
  t <- amt::make_track(i, .x=LON, .y=LAT, .t=DATETIME, crs=sp::CRS("+init=epsg:4326")) 
  r <- track_resample(t, rate=minutes(10), tolerance=minutes(2)) %>% 
    filter_min_n_burst(min_n=3) %>% 
    steps_by_burst() %>% 
    time_of_day(include.crepuscule=FALSE) %>% 
    random_steps(n_control = 100) %>% 
    extract_covariates(substrate, where="both") %>% 
    filter(layer_end != "5") %>% 
    mutate(substrate_start = reclass_substrate(layer_start),
           substrate_end = reclass_substrate(layer_end),
           cos_ta_ = cos(ta_),
           log_sl_ = log(sl_)) %>% 
    filter(!is.infinite(log_sl_)) %>%
    filter(!is.na(ta_))
}
issf034 <- issf(grid_animalshort %>% filter(TRANSMITTER=="BL 034"))


issf051 <- function(i){
  t <- amt::make_track(i, .x=LON, .y=LAT, .t=DATETIME, crs=sp::CRS("+init=epsg:4326")) 
  r <- track_resample(t, rate=minutes(10), tolerance=minutes(2)) %>% 
    filter_min_n_burst(min_n=3) %>% 
    steps_by_burst() %>% 
    time_of_day(include.crepuscule=FALSE) %>% 
    random_steps(n_control = 100) %>% 
    extract_covariates(substrate, where="both") %>% 
    filter(layer_end != "0") %>% 
    filter(layer_end != "0") %>% 
    mutate(substrate_start = reclass_substrate(layer_start),
           substrate_end = reclass_substrate(layer_end),
           cos_ta_ = cos(ta_),
           log_sl_ = log(sl_)) %>% 
    filter(!is.infinite(log_sl_)) %>%
    filter(!is.na(ta_))
}
issf051 <- issf(grid_animalshort %>% filter(TRANSMITTER=="BL 051"))

#---- change substrate reference level ----
issf027 <- mutate(issf027, substrate_start = fct_relevel(substrate_start, "mudL"),
                  substrate_end = fct_relevel(substrate_end, "mudL"))
issf030 <- mutate(issf030, substrate_start = fct_relevel(substrate_start, "mudL"),
                  substrate_end = fct_relevel(substrate_end, "mudL"))
issf032 <- mutate(issf032, substrate_start = fct_relevel(substrate_start, "mudL"),
                  substrate_end = fct_relevel(substrate_end, "mudL"))
issf034 <- mutate(issf034, substrate_start = fct_relevel(substrate_start, "mudL"),
                  substrate_end = fct_relevel(substrate_end, "mudL"))
issf051 <- mutate(issf051, substrate_start = fct_relevel(substrate_start, "mudL"),
                  substrate_end = fct_relevel(substrate_end, "mudL"))

#---- bootstrap ----
set.seed(10)

boot027 <- boot(issf027) 
write.csv(boot027, file = 'boot027.csv')
boot030 <- boot(issf030) 
write.csv(boot030, file = 'boot030.csv')
boot032 <- boot(issf032)
write.csv(boot032, file = 'boot032.csv')
boot034 <- boot(issf034) 
write.csv(boot034, file = 'boot034.csv')
boot051 <- boot(issf051) 
write.csv(boot051, file = 'boot051.csv')

boot027 <- read.csv("boot027.csv")
boot030 <- read.csv("boot030.csv")
boot032 <- read.csv("boot032.csv")
boot034 <- read.csv("boot034.csv")
boot051 <- read.csv("boot051.csv")

#---- RSS plot ----
d027 <- boot027 %>% 
  group_by(term) %>% 
  summarize(
    mean=mean(estimate, na.rm=T),
    ymin=mean-1.96*sd(estimate, na.rm=T),
    ymax=mean+1.96*sd(estimate, na.rm=T)
  )
d027$ID <- "BL 027"
d030 <- boot030 %>% 
  group_by(term) %>% 
  summarize(
    mean=mean(estimate),
    ymin=mean-1.96*sd(estimate),
    ymax=mean+1.96*sd(estimate)
  )
d030$ID <- "BL 030"
d032 <- boot032 %>% 
  group_by(term) %>% 
  summarize(
    mean=mean(estimate),
    ymin=mean-1.96*sd(estimate),
    ymax=mean+1.96*sd(estimate)
  )
d032$ID <- "BL 032"
d034 <- boot034 %>% 
  group_by(term) %>% 
  summarize(
    mean=mean(estimate),
    ymin=mean-1.96*sd(estimate),
    ymax=mean+1.96*sd(estimate)
  )
d034$ID <- "BL 034"
d051 <- boot051 %>% 
  group_by(term) %>% 
  summarize(
    mean=mean(estimate),
    ymin=mean-1.96*sd(estimate),
    ymax=mean+1.96*sd(estimate)
  )
d051$ID <- "BL 051"

#bind together and then do the plot
d <- rbind(d027, d030, d032, d034, d051)
d <- d %>% dplyr::filter(grepl("substrate_end", term))

d$x <- 1:nrow(d)
p1 <- d %>% 
  ggplot(., aes(x = term, y = mean, group = ID, col = ID)) + 
  geom_pointrange(aes(ymin = ymin, ymax = ymax), #indiv points
                  position = position_dodge(width = 0.7), size = 0.8) + 
  geom_hline(yintercept = 0, lty = 2) + #ref class
  labs(x = "Substrate Class", y = "Relative Selection Strength") + 
  theme_light() + 
  scale_x_discrete(labels = c("Mixed sediments with veneer of mud", "No data", "Silt/mud", "Coarse sediments")) +
  coord_flip() 
p1

#---- RSS w bootstrapped df's -----
odds27 <- boot027 %>% split(.$term) %>% purrr::map(~mean(.$estimate, na.rm=T)) %>% bind_rows() %>% 
  mutate(term="mean") %>% gather(key, value, -term)

ssed27 <- exp(-0.11733755) #= the odds are 0.89 times higher that they would choose ssed loc than mudL loc
till27 <- exp(0.43644751) #= the odds are 1.55 times higher that they would choose till loc than mudL loc
csed27 <- exp(0.19389385) #= the odds are 1.21 times higher that they would choose csed loc than mudL loc

odds30 <- boot030 %>% split(.$term) %>% purrr::map(~mean(.$estimate, na.rm=T)) %>% bind_rows() %>% 
  mutate(term="mean") %>% gather(key, value, -term)

ssed30 <- exp(-1.28466419) #= the odds are 0.28 times higher that they would choose ssed loc than mudL loc
till30 <- exp(0.28827541) #= the odds are 1.33 times higher that they would choose till loc than mudL loc
csed30 <- exp(0.03423119) #= the odds are 1.03 times higher that they would choose csed loc than mudL loc


odds32 <- boot032 %>% split(.$term) %>% purrr::map(~mean(.$estimate, na.rm=T)) %>% bind_rows() %>% 
  mutate(term="mean") %>% gather(key, value, -term)

ssed32 <- exp(-0.18337463) #= the odds are 0.83 times higher that they would choose ssed loc than mudL loc
till32 <- exp(0.03465408) #= the odds are 1.04 times higher that they would choose till loc than mudL loc
csed32 <- exp(-0.14353191) #= the odds are 0.87 times higher that they would choose csed loc than mudL loc
nd32 <- exp(0.66813560) #= the odds are 1.95 times higher that they would choose nd loc than mudL loc


odds34 <- boot034 %>% split(.$term) %>% purrr::map(~mean(.$estimate, na.rm=T)) %>% bind_rows() %>% 
  mutate(term="mean") %>% gather(key, value, -term)

ssed34 <- exp(0.4419981) #= the odds are 1.56 times higher that they would choose ssed loc than mudL loc
till34 <- exp(0.2845350) #= the odds are 1.33 times higher that they would choose till loc than mudL loc
csed34 <- exp(0.5205910) #= the odds are 1.68 times higher that they would choose csed loc than mudL loc


odds51 <- boot051 %>% split(.$term) %>% purrr::map(~mean(.$estimate, na.rm=T)) %>% bind_rows() %>% 
  mutate(term="mean") %>% gather(key, value, -term)

csed51 <- exp(0.050629154) #= the odds are 1.05 times higher that they would choose csed loc than mudL loc
till51 <- exp(-0.666625759) #= the odds are 0.51 times higher that they would choose till loc than mudL loc
nd51 <- exp(-0.512658494) #= the odds are 0.60 times higher that they would choose nd loc than mudL loc

#---- effect size w bootstrapping ----
boot027 %>% group_by(term) %>% 
  summarize(q025=quantile(estimate, probs=0.025, na.rm=TRUE), 
            q975=quantile(estimate, probs=0.975, na.rm=TRUE)) %>% 
  mutate(sig=case_when(q025 < 0 & q975 > 0 ~ "no", #is 0 in between q025 and q975?
                       TRUE ~ "yes"))

boot030 %>% group_by(term) %>% 
  summarize(q025=quantile(estimate, probs=0.025, na.rm=TRUE), 
            q975=quantile(estimate, probs=0.975, na.rm=TRUE)) %>% 
  mutate(sig=case_when(q025 < 0 & q975 > 0 ~ "no", #is 0 in between q025 and q975?
                       TRUE ~ "yes"))

boot032 %>% group_by(term) %>% 
  summarize(q025=quantile(estimate, probs=0.025, na.rm=TRUE), 
            q975=quantile(estimate, probs=0.975, na.rm=TRUE)) %>% 
  mutate(sig=case_when(q025 < 0 & q975 > 0 ~ "no", #is 0 in between q025 and q975?
                       TRUE ~ "yes"))

boot034 %>% group_by(term) %>% 
  summarize(q025=quantile(estimate, probs=0.025, na.rm=TRUE), 
            q975=quantile(estimate, probs=0.975, na.rm=TRUE)) %>% 
  mutate(sig=case_when(q025 < 0 & q975 > 0 ~ "no", #is 0 in between q025 and q975?
                       TRUE ~ "yes"))

boot051 %>% group_by(term) %>% 
  summarize(q025=quantile(estimate, probs=0.025, na.rm=TRUE), 
            q975=quantile(estimate, probs=0.975, na.rm=TRUE)) %>% 
  mutate(sig=case_when(q025 < 0 & q975 > 0 ~ "no", #is 0 in between q025 and q975?
                       TRUE ~ "yes"))
