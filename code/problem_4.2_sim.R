library(tidyverse)
library(ivreg)
library(ggeffects)

# functions to run simulation

simulate_data_point = function(alpha1, alpha2, alpha3, beta1, beta2, beta3, Zd_dim, Zs_dim, W_dim) {
  
  # simulate shifters as standard normal
  Zd = rnorm(Zd_dim)
  Zs = rnorm(Zs_dim)
  W = rnorm(W_dim)
  
  # simulate shocks as standard normal
  epsd = rnorm(1)
  epss = rnorm(1)
  
  # solve for log price by imposing supply=demand in equilibrium
  p = (beta2%*%Zs + beta3%*%W + epss - alpha2%*%Zd - alpha3%*%W - epsd)/(alpha1-beta1)
  p = as.numeric(p)
  
  # solve for log quantity
  yd = alpha1%*%p + alpha2%*%Zd + alpha3%*%W + epsd
  yd = as.numeric(yd)
  
  ys = beta1%*%p + beta2%*%Zs + beta3%*%W + epss
  ys = as.numeric(ys)
  
  # check supply = demand at log price p
  check = as.logical(abs(yd-ys) < 1e-6)
  
  # return (price, quantity, (observables))
  if (check){
    return(list(
      "price" = p, 
      "quantity" = yd,
      list("demand_shifter" = Zd, "supply_shifter" = Zs, "common_shifter" = W)
      ))
  } else {
    print("Could not compute market-clearing price")
    return(NULL)
  }
}

simulate_data = function(N, demand_coefs, supply_coefs) {
  
  # unpack arguments
  alpha1 = demand_coefs[[1]]
  alpha2 = demand_coefs[[2]]
  alpha3 = demand_coefs[[3]]
  
  beta1 = supply_coefs[[1]]
  beta2 = supply_coefs[[2]]
  beta3 = supply_coefs[[3]]
  
  # dimensions
  Zd_dim = length(alpha2)
  Zs_dim = length(beta2)
  W_dim = length(beta3)
  
  # prepopulate matrices
  Y_vec = rep(NA, N)
  p_vec = rep(NA, N)
  Zd_mat = matrix(NA, nrow = N, ncol = Zd_dim)
  Zs_mat = matrix(NA, nrow = N, ncol = Zs_dim)
  W_mat =  matrix(NA, nrow = N, ncol = W_dim)
  
  # run simulations
  for (i in 1:N) {
    data_point_sim = simulate_data_point(alpha1, alpha2, alpha3, beta1, beta2, beta3, Zd_dim, Zs_dim, W_dim)
    
    # fill out matrices
    Y_vec[i] = data_point_sim$quantity[1]
    p_vec[i] = data_point_sim$price[1]
    Zd_mat[i,] = data_point_sim[[3]]$demand_shifter
    Zs_mat[i,] = data_point_sim[[3]]$supply_shifter
    W_mat[i,] = data_point_sim[[3]]$common_shifter
  
  }
  
  return(list(
    "quantities" = Y_vec,
    "prices" = p_vec,
    "demand_shifters" = Zd_mat,
    "supply_shifters" = Zs_mat,
    "common_shifters" = W_mat
  ))
}

# simulate a dataset
set.seed(14382)
N = 10^3
demand_coefficients = list(
  "alpha1" = -0.8, # demand decreasing in price
  "alpha2" = c(1,.5), # can have multidimensional shifters,
  "alpha3" = c(.5)
)
supply_coefficients = list(
  "beta1" = 1.1, # supply increasing in price
  "beta2" = c(1.4,0.3,0.1),
  "beta3" = c(.3)
)

data_list = simulate_data(N, demand_coefficients, supply_coefficients)

## plot prices and quantities

# ols
ols_df = tibble(
  "p" = data_list$prices, 
  "y" = data_list$quantities,
  "w" = data_list$common_shifters[,1]
)

ols_model = lm(
  formula = y ~ p + w,
  data = ols_df
)

ols_elasticity = summary(ols_model)$coefficients[2,"Estimate"]
ols_se = sqrt(diag(vcovHC(ols_model, type = "HC1")))[2]

summary(ols_model, diagnostics=TRUE)

ols_pred = ggpredict(ols_model, terms = "p")
ols_pred$model = "OLS"


# 2sls demand
demand_df = tibble(
  "p" = data_list$prices, 
  "y" = data_list$quantities,
  "w" = data_list$common_shifters,
  "z1" = data_list$supply_shifters[,1],
  "z2" = data_list$supply_shifters[,2],
  "z3" = data_list$supply_shifters[,3]
)
  
demand_2sls_model = ivreg(
  formula = y ~ p + w | z1 +z2 + z3,
  data = demand_df
)

demand_2sls_elasticity = summary(demand_2sls_model)$coefficients[2,"Estimate"]
demand_2sls_se = sqrt(diag(vcovHC(demand_2sls_model, type = "HC1")))[2]

summary(demand_2sls_model, diagnostics=TRUE)

demand_pred = ggpredict(demand_2sls_model, terms = "p")
demand_pred$model = "Demand Curve 2SLS"

# 2sls supply
supply_df = tibble(
  "p" = data_list$prices, 
  "y" = data_list$quantities,
  "w" = data_list$common_shifters,
  "z1" = data_list$demand_shifters[,1],
  "z2" = data_list$demand_shifters[,2]
)

supply_2sls_model = ivreg(
  formula = y ~ p + w | z1 + z2,
  data = supply_df
)

summary(supply_2sls_model, diagnostics=TRUE)

supply_2sls_elasticity = summary(supply_2sls_model)$coefficients[2,"Estimate"]
supply_2sls_se = sqrt(diag(vcovHC(supply_2sls_model, type = "HC1")))[2]

supply_pred = ggpredict(supply_2sls_model, terms = "p")
supply_pred$model = "Supply Curve 2SLS"

# scatter

pred_all <- rbind(ols_pred, demand_pred, supply_pred)

ggplot() +
  geom_point(data = tibble("price" = data_list$prices, "quantity" = data_list$quantities), aes(x = price, y = quantity), size = 0.6) +
  geom_line(data = pred_all, aes(x=x, y = predicted, color = model)) +
  geom_ribbon(data = pred_all, aes(x=x, ymin = conf.low, ymax = conf.high, fill = model), alpha = 0.2) +
  theme_minimal() + 
  xlim(-5,5) +
  ylim(-5,5) +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    title = "Supply and Demand System Simulation",
    x = "log Price",
    y = "log Quantity"
  )

ggsave(
  filename = "output/problem_4_2_plot.png",    # file path + name
  plot = last_plot(),          # plot object
  device = "png",              # format: "png", "pdf", "jpeg", "tiff", "svg"
  width = 8,                   # width in inches
  height = 5,                  # height in inches
  units = "in",                # units: "in", "cm", "mm"
  dpi = 300,                   # resolution: 300dpi is standard for print
  bg = "white"                 # background color, important for PNG
)

