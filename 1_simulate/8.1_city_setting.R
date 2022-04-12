## populations
## population of susceptible

df_raw <- read.xlsx('data/demo.xlsx', sheet = city_list[c])
df_raw$total <- df_raw$total * (10^4)
pop_s <- df_raw |> 
  filter(value != 'total') |> 
  mutate(full_vaccine = total * (full_vaccine_cover - booster_vaccine_cover),
         booster_vaccine = total * booster_vaccine_cover,
         un_vaccine = total - full_vaccine - booster_vaccine) |> 
  select(un_vaccine, full_vaccine, booster_vaccine) |> 
  as.matrix() |> 
  matrix(nrow = 1) |> 
  round()
pop_n <- pop_s
pop_e <- matrix(rep(0, length(pop_s)), nrow = 1)
remove(df_raw)

### population of pre-symptomatic
df_inf <- read.xlsx('data/infections income.xlsx', sheet = 'presymptomatic')
pop_ip <- df_inf |> 
  select(-value) |> 
  as.matrix() |> 
  matrix(nrow = 1)
### population of symptomatic
df_inf <- read.xlsx('data/infections income.xlsx', sheet = 'symptomatic')
pop_is <- df_inf |> 
  select(-value) |> 
  as.matrix() |> 
  matrix(nrow = 1)
### population of asymptomatic
df_inf <- read.xlsx('data/infections income.xlsx', sheet = 'asymptomatic')
pop_ia <- df_inf |> 
  select(-value) |> 
  as.matrix() |> 
  matrix(nrow = 1)
remove(df_inf)

### population of remove/recovery
pop_r <- matrix(rep(0, length(pop_s)), nrow = 1)

## read contact matrix

df_cm <- read.xlsx('data/contact_before.xlsx', sheet = city_list[c])
contact_mt_before <- df_cm |> select(-value) |> as.matrix()
df_cm <- read.xlsx('data/contact_after.xlsx', sheet = city_list[c])
contact_mt_after <- df_cm |>   select(-value) |> as.matrix()
remove(df_cm)

mcmc_para <- data.frame(
  omega_1 = 1/rgamma(n, shape = df_para[1, 'value_shape'], rate = df_para[1, 'value_rate']),
  omega_2 = 1/rgamma(n, shape = df_para[2, 'value_shape'], rate = df_para[2, 'value_rate']),
  omega_3 = 1/rgamma(n, shape = df_para[3, 'value_shape'], rate = df_para[3, 'value_rate'])
  # omega_1 = 1/(df_para[1, 'value_shape']/ df_para[1, 'value_rate']),
  # omega_2 = 1/(df_para[2, 'value_shape']/ df_para[2, 'value_rate']),
  # omega_3 = 1/(df_para[3, 'value_shape']/ df_para[3, 'value_rate'])
  # omega_3 = 1/1
) |> 
  mutate(constant = ((1-percent_asym)*omega_2 + percent_asym*omega_1)/
           ((1-percent_asym)*omega_2*rb_is/gamma_2 + 
              percent_asym*omega_1*rb_ia/gamma_1 + 
              (1-percent_asym)*omega_2*rb_ip/omega_3),
         beta6 = compute.beta(RP_numbers[1], constant, contact_mt_before),
         beta8 = compute.beta(RP_numbers[2], constant, contact_mt_before),
         beta10 = compute.beta(RP_numbers[3], constant, contact_mt_before),
         beta12 = compute.beta(RP_numbers[4], constant, contact_mt_before)) |> 
  pivot_longer(cols = -c(omega_1, omega_2, omega_3, constant), values_to = 'beta', names_to = 'RP_number') |> 
  mutate(RP_number = as.numeric(gsub('beta', '', RP_number))) |> 
  select(omega_1, omega_2, omega_3, beta, RP_number)
