library(shiny)

options(shiny.maxRequestSize = 50 * 1024^2)

# ==============================================================================
# 1. MATEMATIKA A GENERÁTORY (BEZ ZMIEN)
# ==============================================================================

generuj_logseriove <- function(S, theta) {
  k <- 1:S; vahy <- (theta^k)/k
  pocty <- sample(k, size=S, replace=TRUE, prob=vahy)
  return(rep(1:length(pocty), pocty))
}

generuj_geometricke <- function(S, p) {
  k <- 1:S
  vahy <- p * (1 - p)^(k - 1) 
  pocty <- sample(k, size = S, replace = TRUE, prob = vahy)
  return(rep(1:length(pocty), pocty))
}

generuj_lognormalne <- function(S, meanlog, sdlog) {
  abundancie <- rlnorm(n = S, meanlog = meanlog, sdlog = sdlog)
  abundancie <- ceiling(abundancie)
  return(rep(1:S, abundancie))
}

generuj_macarthur_stick <- function(S, N_total = NULL) {
  if (is.null(N_total)) N_total <- S * 20
  body_zlomu <- sort(runif(S - 1))
  vsetky_body <- c(0, body_zlomu, 1)
  podiely <- diff(vsetky_body)
  abundancie <- round(podiely * N_total)
  abundancie[abundancie == 0] <- 1
  return(rep(1:S, abundancie))
}

generuj_yule <- function(S, rho) {
  k <- 1:S 
  vahy <- rho * beta(k, rho + 1)
  pocty <- sample(k, size = S, replace = TRUE, prob = vahy)
  return(rep(1:length(pocty), pocty))
}

vypocitaj_alfa <- function(S, N) {
  if (S <= 0 || N <= S) return(NA)
  rovnica <- function(a) S - a * log(1 + N/a)
  tryCatch({ uniroot(rovnica, interval=c(1e-5, 1e10))$root }, error=function(e) NA)
}

vypocitaj_fisher_parametre <- function(vzorka) {
  S <- length(unique(vzorka))
  N <- length(vzorka)
  if (S <= 0 || N <= S) return(list(alpha = NA, k = 0, use_k = FALSE))
  freq <- as.numeric(table(vzorka))
  if (length(freq) < 2) {
    alfa <- vypocitaj_alfa(S, N)
    return(list(alpha = alfa, k = 0, use_k = FALSE))
  }
  a_n <- table(freq)
  n_vals <- as.numeric(names(a_n)); a_vals <- as.numeric(a_n)
  term1 <- sum(n_vals * (n_vals - 1) * a_vals)
  alfa <- vypocitaj_alfa(S, N)
  if (is.na(alfa)) return(list(alpha = NA, k = 0, use_k = FALSE))
  x <- N / (N + alfa); term2 <- (S^2) / (2 * alfa)
  discrepancy <- term1 - term2; log_N_over_S <- log10(N/S)
  log_vals <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5)
  i_vals <- c(0.1971, 0.2882, 0.3914, 0.5054, 0.6295, 0.7639, 0.9076, 1.0608, 1.2232, 1.3950,
              1.5762, 1.7665, 1.9661, 2.1751, 2.3934, 2.6211, 2.8582, 3.1047, 3.3606, 3.6260,
              3.9009, 4.1854, 4.4791, 4.7825, 5.0954, 5.4178, 5.7498, 6.0912, 6.4421, 6.8026, 7.1726, 7.5521)
  if (log_N_over_S < 0.4) { i_over_S <- 0.1 } else if (log_N_over_S > 3.5) { i_over_S <- 8.0 } else {
    i_over_S <- approx(log_vals, i_vals, xout = log_N_over_S, rule = 2)$y
  }
  i <- i_over_S * S; se <- sqrt(i)
  if (!is.na(se) && se > 0 && abs(discrepancy) > 2 * se) {
    k_estimate <- discrepancy / i
    if (k_estimate > 0 && k_estimate < 10) {
      alpha_k <- S * k_estimate
      return(list(alpha = alpha_k, k = k_estimate, use_k = TRUE))
    }
  }
  return(list(alpha = alfa, k = 0, use_k = FALSE))
}

predikuj_fisher <- function(vzorka, cas1, cas2) {
  params <- vypocitaj_fisher_parametre(vzorka)
  if (is.na(params$alpha)) return(NA)
  S <- length(unique(vzorka)); N <- length(vzorka); N_new <- N * (cas2 / cas1)
  if (params$use_k && params$k > 0) {
    alfa <- params$alpha; k <- params$k
    odhad <- (alfa / k) * ((alfa / (N + alfa))^k - (alfa / (N + N_new + alfa))^k)
    if (is.na(odhad) || is.infinite(odhad) || odhad < 0) {
      alfa <- vypocitaj_alfa(S, N); if (is.na(alfa)) return(NA)
      return(alfa * log(1 + N_new / (N + alfa)))
    }
    return(odhad)
  } else {
    alfa <- params$alpha; return(alfa * log(1 + N_new / (N + alfa)))
  }
}

predikuj_smoothed_gt <- function(vzorka, cas1, cas2, k_param = 20) {
  t <- cas2 / cas1; freq <- as.numeric(table(vzorka)); max_i <- max(freq)
  phi <- numeric(max_i); for(i in 1:max_i) phi[i] <- sum(freq == i)
  q <- 1 / (1 + t); odhad <- 0; limit <- min(length(phi), k_param)
  for(i in 1:limit) {
    if (phi[i] > 0) {
      prob_L_geq_i <- pbinom(i - 1, k_param, q, lower.tail = FALSE)
      odhad <- odhad + ((-1)^(i + 1) * (t^i) * prob_L_geq_i * phi[i])
    }
  }
  return(odhad)
}

predikuj_gt <- function(vzorka, cas1, cas2) {
  t <- cas2 / cas1
  if (t > 1) return(predikuj_smoothed_gt(vzorka, cas1, cas2))
  freq <- as.numeric(table(vzorka)); max_i <- max(freq); phi <- numeric(max_i)
  for(i in 1:max_i) phi[i] <- sum(freq == i)
  odhad <- 0; limit <- min(length(phi), 50)
  for(i in 1:limit) { if(phi[i] > 0) odhad <- odhad + (-((-t)^i) * phi[i]) }
  return(odhad)
}

predikuj_subsampling_weighted <- function(vzorka, cas1, cas2, n_iter = 50) {
  n1 <- length(vzorka); if (n1 < 2) return(NA)
  vaha <- cas1 / (cas1 + cas2)
  n_trening <- round(vaha * n1); n_test <- round((1 - vaha) * n1)
  if (n_trening < 1) n_trening <- 1; if (n_test < 1) n_test <- 1; if (n_trening >= n1) n_trening <- n1 - 1
  odhady <- numeric(n_iter)
  for (i in 1:n_iter) {
    idx <- sample(1:n1, n_trening, replace = FALSE)
    sub_trening <- vzorka[idx]; sub_test <- vzorka[-idx]
    if(length(sub_test) > n_test) sub_test <- sub_test[1:n_test]
    odhady[i] <- length(setdiff(unique(sub_test), unique(sub_trening)))
  }
  return(mean(odhady))
}


predikuj_subsampling_analyticky <- function(vzorka, cas1, cas2) {
  n <- length(vzorka)
  if (n < 2) return(NA)
  m <- round(n * (cas2 / cas1))
  if (m < 1) m <- 1
  
  vaha <- cas1 / (cas1 + cas2)
  n1 <- round(vaha * n)
  if (n1 < 1) n1 <- 1
  if (n1 >= n) n1 <- n - 1
  n2 <- n - n1
  
  freq <- table(table(vzorka))
  j_vals <- as.numeric(names(freq))
  a_j <- as.numeric(freq)
  
  # ---- Krok 1: výpočet hat U_w cez log-binomické koeficienty ----
  U_sub <- 0
  for (i in 1:length(j_vals)) {
    j <- j_vals[i]
    if (j <= n2) {
      p <- exp(lchoose(n - j, n1) - lchoose(n, n1))
      U_sub <- U_sub + a_j[i] * p
    }
  }
  
  if (U_sub <= 0) return(0)
  
  # ---- Krok 2: pomer spomaľovania r ----
  a1 <- sum(a_j[j_vals == 1])
  sklon_koniec <- a1 / n
  sklon_priemer <- U_sub / n2
  
  if (sklon_priemer <= 0) return(U_sub)
  r <- sklon_koniec / sklon_priemer
  
  # ---- Krok 3: numericky stabilná teoretická funkcia r(β) ----
  # Použijeme parametrizáciu cez ρ = n1/n; pre β > 50 prejdeme
  # na asymptotický log-výpočet, aby sme sa vyhli pretečeniu.
  rho <- n1 / n
  
  teoreticky_r <- function(beta) {
    if (abs(1 - beta) < 1e-10) {
      return((1 - rho) / (-log(rho)))
    }
    if (beta > 50) {
      # asymptotika: r(β) ≈ (1-ρ)(β-1)·ρ^(β-1)
      log_r <- log(1 - rho) + log(beta - 1) + (beta - 1) * log(rho)
      return(if (log_r > -700) exp(log_r) else 0)
    }
    return((1 - rho) * (1 - beta) / (1 - rho^(1 - beta)))
  }
  
  # ---- Krok 4: ADAPTÍVNE HRANICE BISEKCIE ----
  # β_lower zodpovedá r = 0.999 (takmer žiadne spomaľovanie)
  # β_upper zodpovedá r = 0.001 (takmer úplná saturácia)
  # Vonkajší bracket [1e-8, 1000] je dostatočne široký pre všetky ρ
  # a r(β) je striktne klesajúca, takže koreň vždy existuje.
  
  beta_lower <- tryCatch({
    uniroot(function(b) teoreticky_r(b) - 0.999,
            interval = c(1e-8, 1000), tol = 1e-8)$root
  }, error = function(e) 1e-4)
  
  beta_upper <- tryCatch({
    uniroot(function(b) teoreticky_r(b) - 0.001,
            interval = c(1e-8, 1000), tol = 1e-8)$root
  }, error = function(e) 50)
  
  # ---- Krok 5: hlavná bisekcia s adaptívnymi hranicami ----
  # Ošetrenie prípadov, keď pozorované r leží mimo [0.001, 0.999]:
  if (r >= teoreticky_r(beta_lower)) {
    beta_odhad <- beta_lower
  } else if (r <= teoreticky_r(beta_upper)) {
    beta_odhad <- beta_upper
  } else {
    beta_odhad <- tryCatch({
      uniroot(function(b) teoreticky_r(b) - r,
              interval = c(beta_lower, beta_upper),
              tol = 1e-6)$root
    }, error = function(e) 1)
  }
  
  # ---- Krok 6: škálovací faktor R(β̂) ----
  if (abs(1 - beta_odhad) < 1e-10) {
    R_faktor <- log((n + m) / n) / log(n / n1)
  } else {
    b <- 1 - beta_odhad
    R_faktor <- ((n + m)^b - n^b) / (n^b - n1^b)
  }
  
  if (is.na(R_faktor) || R_faktor < 0) R_faktor <- 1
  
  return(U_sub * R_faktor)
}

zisti_najlepsi_model_aic <- function(vzorka) {
  x <- as.numeric(table(vzorka)); N <- sum(x); S <- length(x)
  
  # Log-Séria
  alfa <- vypocitaj_alfa(S, N); if(is.na(alfa)) alfa <- S
  p_ls <- N / (N + alfa); const_ls <- -1 / log(1 - p_ls)
  ll_ls <- sum(log(const_ls * (p_ls^x) / x + 1e-10))
  aic_ls <- 4 - 2*ll_ls  
  
  # Geometrické
  p_geom <- S/N; probs_geom <- dgeom(x - 1, prob = p_geom)
  ll_geom <- sum(log(probs_geom + 1e-10)); aic_geom <- 2 - 2*ll_geom 
  
  # Log-Normálne
  log_x <- log(x); mu <- mean(log_x); sigma <- sd(log_x)
  probs_ln <- dlnorm(x, meanlog = mu, sdlog = sigma)
  ll_ln <- sum(log(probs_ln + 1e-10)); aic_ln <- 4 - 2*ll_ln  
  
  
  
  # Yule - vahy podľa rho * beta(k, rho+1)
  ll_yl <- tryCatch({
    # Pre Yule rozdelenie odhadneme rho metódou momentov: priemerná hodnota ~ rho/(rho-1)
    # Použijeme jednoduchú MLE aproximáciu
    rang <- rank(-x, ties.method = "first")
    rho_seq <- seq(0.5, 5, by = 0.1)
    best_ll <- -Inf
    for (rho in rho_seq) {
      vahy <- rho * beta(rang, rho + 1)
      vahy <- vahy / sum(vahy)
      vahy[vahy <= 0] <- 1e-10
      ll_test <- sum(x * log(vahy + 1e-10))
      if (ll_test > best_ll) best_ll <- ll_test
    }
    best_ll
  }, error = function(e) -Inf)
  aic_yl <- 2 - 2*ll_yl
  
  vysledky <- c("Log-Séria"=aic_ls, "Geometrické"=aic_geom, "Log-Normálne"=aic_ln, 
                "Yule"=aic_yl)
  return(names(which.min(vysledky)))
}

# --- UPRAVENÁ FUNKCIA S ITERÁCIAMI ---
analyzuj_jadro_simulacie <- function(populacia, cas_tr, cas_pr, model, iteracie = 10) {
  N_tot <- length(populacia); rate <- 0.1 
  n1 <- round(N_tot * rate * cas_tr); n2 <- round(N_tot * rate * cas_pr)
  if(n1+n2 > N_tot) { sc<-N_tot/(n1+n2); n1<-floor(n1*sc); n2<-floor(n2*sc) }
  
  # Príprava vektorov pre výsledky z každej iterácie
  res_pf <- numeric(iteracie); res_pg <- numeric(iteracie)
  res_ps <- numeric(iteracie); res_pa <- numeric(iteracie)
  res_real <- numeric(iteracie)
  
  # Beh iterácií
  for(i in 1:iteracie) {
    idx <- sample(1:N_tot)
    v1 <- populacia[idx[1:n1]]
    v2 <- populacia[idx[(n1+1):(n1+n2)]]   
    
    res_real[i] <- length(setdiff(unique(v2), unique(v1)))
    res_pf[i] <- predikuj_fisher(v1, cas_tr, cas_pr)
    res_pg[i] <- predikuj_gt(v1, cas_tr, cas_pr)
    res_ps[i] <- predikuj_subsampling_weighted(v1, cas_tr, cas_pr)
    res_pa[i] <- predikuj_subsampling_analyticky(v1, cas_tr, cas_pr)
  }
  
  # Spriemerovanie
  mean_real <- mean(res_real, na.rm = TRUE)
  mean_pf <- mean(res_pf, na.rm = TRUE)
  mean_pg <- mean(res_pg, na.rm = TRUE)
  mean_ps <- mean(res_ps, na.rm = TRUE)
  mean_pa <- mean(res_pa, na.rm = TRUE)
  
  # Bezpečný výpočet chyby z priemerov
  err <- function(p,r) if(is.na(p) || is.na(r) || r==0) 0 else (p-r)/r*100
  
  data.frame(Model = model, 
             Pred_Fisher = round(mean_pf, 1), 
             Pred_GT = round(mean_pg, 1), 
             Pred_Sub = round(mean_ps, 1), 
             Pred_Sub_Analyticky = round(mean_pa, 1), 
             REALITA = round(mean_real, 1), 
             Chyba_Fisher = round(err(mean_pf, mean_real), 1), 
             Chyba_GT = round(err(mean_pg, mean_real), 1), 
             Chyba_Sub = round(err(mean_ps, mean_real), 1), 
             Chyba_Sub_Analyticky = round(err(mean_pa, mean_real), 1))
}

# ==============================================================================
# 2. DIZAJN APLIKÁCIE 
# ==============================================================================

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { font-family: 'Segoe UI', sans-serif; background-color: #f4f7f6; }
      .centered-container { display: flex; flex-direction: column; align-items: center; justify-content: center; height: 90vh; text-align: center; }
      .btn-custom { border-radius: 30px; padding: 20px 40px; font-size: 20px; font-weight: bold; border: none; box-shadow: 0 8px 15px rgba(0, 0, 0, 0.1); transition: all 0.3s ease 0s; margin: 20px; color: white; }
      .btn-custom:hover { transform: translateY(-5px); box-shadow: 0 15px 20px rgba(0, 0, 0, 0.2); }
      .btn-science { background: linear-gradient(45deg, #0e6332, #026112); }
      .panel-science { border-top: 5px solid #0e6332; background: white; padding: 20px; border-radius: 10px; box-shadow: 0 4px 10px rgba(0,0,0,0.05); }
      .btn-sim { background: linear-gradient(45deg, #3498db, #2980b9); }
      .panel-sim { border-top: 5px solid #3498db; background: white; padding: 20px; border-radius: 10px; box-shadow: 0 4px 10px rgba(0,0,0,0.05); }
      .result-box { background: white; padding: 20px; border-radius: 8px; margin-top: 20px; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }
      .btn-back { margin-top: 20px; background-color: #95a5a6; color: white; border-radius: 20px; }
      .input-centered { max-width: 600px; margin: 0 auto; text-align: left; }
      .btn-remove-file { background-color: #e74c3c; color: white; border: none; font-size: 12px; padding: 5px 10px; border-radius: 4px; margin-bottom: 15px; cursor: pointer;}
      .btn-remove-file:hover { background-color: #c0392b; }
      .table-container { overflow-x: auto; width: 100%; }
      
      /* Štýl pre download box */
      .download-box {
        margin-top: 20px;
        padding: 15px;
        background: #f8fbfd;
        border-radius: 8px;
        border: 1px solid #e1e8ed;
      }
      .download-box h4 {
        margin-top: 0;
        margin-bottom: 12px;
        
      }
      
      /* Tmavomodré tlačidlo */
      .btn-darkblue {
        background: linear-gradient(45deg, #0a1f44, #001a5e) !important;
        color: white !important;
        border: none !important;
        font-weight: bold;
        width: 100%;
        margin-bottom: 8px;
        padding: 6px 12px;
        border-radius: 4px;
        text-align: center;
        display: inline-block;
        text-decoration: none;
      }
      .btn-darkblue:hover {
        background: linear-gradient(45deg, #001340, #000c3d) !important;
        color: white !important;
        text-decoration: none;
      }
      .btn-darkblue:focus, .btn-darkblue:active {
        color: white !important;
        text-decoration: none;
      }
      
      
      /* Štýl pre tlačidlo Rozšírené nastavenia */
      .btn-advanced-toggle {
        background: none;
        border: none;
        color: #7f8c8d;
        font-size: 14px;
        text-align: center;
        width: 100%;
        padding: 8px;
        cursor: pointer;
        text-decoration: underline;
        transition: color 0.2s ease;
      }
      .btn-advanced-toggle:hover {
        color: #2980b9;
      }
      .btn-advanced-toggle .arrow {
        display: inline-block;
        transition: transform 0.4s ease;
        margin-left: 5px;
      }
      .btn-advanced-toggle.open .arrow {
        transform: rotate(180deg);
      }
      
      /* Animácia rozbaľovania */
      .advanced-panel {
        max-height: 0;
        overflow: hidden;
        transition: max-height 0.5s ease-in-out, padding 0.5s ease-in-out, margin 0.5s ease-in-out, opacity 0.4s ease-in-out;
        opacity: 0;
        padding: 0 15px;
        margin-top: 0;
        background: #f8fbfd;
        border-radius: 8px;
        border: 1px solid #e1e8ed;
      }
      .advanced-panel.open {
        max-height: 1500px;
        opacity: 1;
        padding: 15px;
        margin-top: 10px;
      }
      .advanced-panel h4 {
        color: #2980b9;
        margin-top: 10px;
        margin-bottom: 8px;
        border-bottom: 1px solid #d6e4ee;
        padding-bottom: 5px;
      }
    ")),
    tags$script(HTML("
      $(document).on('click', '#toggle_advanced', function() {
        $('#advanced_panel').toggleClass('open');
        $(this).toggleClass('open');
        var isOpen = $(this).hasClass('open');
        $(this).find('.toggle-text').text(isOpen ? 'Skryť rozšírené nastavenia' : 'Rozšírené nastavenia');
      });
    "))
  ),
  tabsetPanel(id = "main_nav", type = "hidden",
              tabPanel("home",
                       div(class = "centered-container",
                           h1("Problem odhadu počtu nepozorovaných druhov", style = "margin-bottom: 50px; color: #2c3e50;"),
                           div(actionButton("go_science", "Vedecká Predpoveď", class = "btn-custom btn-science"),
                               actionButton("go_sim", "Simulácia Metód", class = "btn-custom btn-sim"))
                       )
              ),
              tabPanel("science",
                       div(style = "padding: 20px;",
                           div(class = "panel-science input-centered",
                               h2("Vstup dát", style="text-align: center; color: #026112;"),
                               uiOutput("file_upload_ui"),
                               actionButton("remove_file", "X Odstrániť súbor", class = "btn-remove-file"),
                               textAreaInput("manual_sc", "Alebo manuálne (oddelené čiarkou)", placeholder = "Vlk, Vlk, Rys, 10, 20...", rows=2),
                               fluidRow(column(6, numericInput("t1_sc", "Čas zberu (t1)", value = 12, min = 0.001)),
                                        column(6, numericInput("t2_sc", "Budúci čas (t2)", value = 6, min = 0.001))),
                               actionButton("run_sc", "Analyzovať dáta", class = "btn-custom btn-science", style = "width: 100%; margin: 10px 0;")
                           ),
                           uiOutput("results_sc"),
                           actionButton("back_home_1", "Späť na úvod", class = "btn btn-back")
                       )
              ),
              tabPanel("simulation",
                       div(style = "padding: 20px;",
                           div(class = "panel-sim input-centered",
                               h2("Parametre Simulácie", style="text-align: center; color: #2980b9;"),
                               fluidRow(
                                 column(6, numericInput("S_sim", "Počet druhov (S)", value = 1000, min = 100)),
                                 column(6, numericInput("iter_sim", "Počet iterácií", value = 10, min = 1))
                               ),
                               fluidRow(
                                 column(6, numericInput("t1_sim", "Tréning (t1)", value = 1, min = 0.1)),
                                 column(6, numericInput("t2_sim", "Test (t2)", value = 1, min = 0.1))
                               ),
                               actionButton("run_sim", "Spustiť Simuláciu", class = "btn-custom btn-sim", style = "width: 100%; margin: 10px 0;"),
                               
                               # --- TLAČIDLO ROZŠÍRENÉ NASTAVENIA ---
                               tags$button(id = "toggle_advanced", class = "btn-advanced-toggle",
                                           tags$span(class = "toggle-text", "Rozšírené nastavenia"),
                                           tags$span(class = "arrow", HTML("▼"))
                               ),
                               
                               # --- ROZBAĽOVACÍ PANEL S PARAMETRAMI DISTRIBÚCIÍ ---
                               div(id = "advanced_panel", class = "advanced-panel",
                                   h4("Log-Séria"),
                                   numericInput("param_ls_theta", "θ (theta)", value = 0.95, min = 0.01, max = 0.999, step = 0.01),
                                   
                                   h4("Log-Normálne"),
                                   fluidRow(
                                     column(6, numericInput("param_ln_meanlog", "meanlog (μ)", value = 3.0, step = 0.1)),
                                     column(6, numericInput("param_ln_sdlog", "sdlog (σ)", value = 1.0, min = 0.01, step = 0.1))
                                   ),
                                   
                                   h4("Geometrické"),
                                   numericInput("param_gm_p", "p (pravdepodobnosť)", value = 0.1, min = 0.001, max = 0.999, step = 0.01),
                                   
                                   h4("MacArthur (Broken Stick)"),
                                   numericInput("param_bs_N", "N_total (celkový počet jedincov)", value = 40000, min = 100, step = 1000),
                                   
                                   h4("Yule"),
                                   numericInput("param_yl_rho", "ρ (rho)", value = 1.5, min = 0.1, step = 0.1)
                               )
                           ),
                           uiOutput("results_sim"),
                           actionButton("back_home_2", "Späť na úvod", class = "btn btn-back")
                       )
              )
  )
)

# ==============================================================================
# 3. SERVER
# ==============================================================================

server <- function(input, output, session) {
  
  # --- POMOCNÉ FUNKCIE NA VYKRESLENIE GRAFOV (na opätovné použitie pri sťahovaní) ---
  vykresli_model_sim <- function(p, t, t1, t2, rate = 0.1) {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    cnt <- sort(table(p), decreasing = TRUE)
    sad <- table(cnt)
    x <- as.numeric(names(sad)); y <- as.numeric(sad)
    df <- data.frame(x, y); df <- df[order(df$x), ]; df <- df[df$x <= 20, ]
    if (nrow(df) == 0) df <- data.frame(x = 1, y = 0)
    barplot(df$y, names.arg = df$x, col = "#3498db", border = NA,
            main = paste("Distribúcia abundancií modelu", t),
            xlab = "Počet jedincov (n)", ylab = "Počet druhov (S_n)")
    
    N_tot <- length(p)
    n1 <- round(N_tot * rate * t1); n2 <- round(N_tot * rate * t2)
    if (n1 + n2 > N_tot) { scale <- N_tot / (n1 + n2); n1 <- floor(n1 * scale); n2 <- floor(n2 * scale) }
    cely_lov <- sample(p)
    je_novy <- !duplicated(cely_lov)
    kumulativny_pocet <- cumsum(je_novy)
    plot(1:length(kumulativny_pocet), kumulativny_pocet, type = "l", col = "gray", lwd = 1,
         main = paste("Krivka objavovania modelu", t),
         xlab = "Počet jedincov", ylab = "Kumulatívny počet druhov")
    if (n1 > 0) lines(1:n1, kumulativny_pocet[1:n1], col = "#0dc4fc", lwd = 2)
    if (n2 > 0 && (n1 + n2) <= length(kumulativny_pocet)) {
      lines((n1 + 1):(n1 + n2), kumulativny_pocet[(n1 + 1):(n1 + n2)], col = "#e67e22", lwd = 2)
    }
    abline(v = n1, lty = 2, col = "black")
    legend("bottomright", legend = c("Tréning", "Test", "Zvyšok"),
           fill = c("#0dc4fc", "#e67e22", "gray"), cex = 0.6, bty = "n")
  }
  
  vykresli_graf_sc <- function(v) {
    if (is.null(v) || length(v) == 0) return(NULL)
    freq <- table(v); counts <- sort(as.numeric(freq), decreasing=T); labels_x <- names(sort(freq, decreasing=TRUE))
    if(length(counts) == 0) return(NULL)
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
    top_n <- min(15, length(counts)); barplot(counts[1:top_n], names.arg=labels_x[1:top_n], las=2, col="#0e6332", border=NA, main="Dominancia (Top 15)", ylab="Jedincov", cex.names=0.8)
    sad <- table(counts); x<-as.numeric(names(sad)); y<-as.numeric(sad); df<-data.frame(x,y); df<-df[df$x<=30,]; if(nrow(df)==0) df<-data.frame(x=1,y=0)
    barplot(df$y, names.arg=df$x, col="#026112", border=NA, main="Distribúcia druhovej abundancie", xlab="Jedincov", ylab="Druhov")
    je_novy <- !duplicated(v); kumulativny_pocet <- cumsum(je_novy)
    plot(1:length(v), kumulativny_pocet, type="l", lwd=2, col="#27ae60", main="Krivka objavovania druhov", xlab="Počet pozorovaní", ylab="Kumulatívny počet druhov")
    grid(col="lightgray", lty=1); points(length(v), kumulativny_pocet[length(v)], col="red", pch=19, cex=1.5)
  }
  
  reset_trigger <- reactiveVal(0)
  rv <- reactiveValues(data = NULL, analyzed = FALSE, t1 = 12, t2 = 6)
  
  observeEvent(input$go_science, { updateTabsetPanel(session, "main_nav", "science") })
  observeEvent(input$go_sim, { updateTabsetPanel(session, "main_nav", "simulation") })
  observeEvent(input$back_home_1, { updateTabsetPanel(session, "main_nav", "home") })
  observeEvent(input$back_home_2, { updateTabsetPanel(session, "main_nav", "home") })
  
  output$file_upload_ui <- renderUI({
    reset_trigger()
    fileInput("file_sc", "Nahrať súbor (.csv, .txt)", accept = c(".csv", ".txt"))
  })
  
  observeEvent(input$remove_file, {
    reset_trigger(reset_trigger() + 1)
    rv$data <- NULL
    rv$analyzed <- FALSE
    updateTextAreaInput(session, "manual_sc", value = "")
  })
  
  observeEvent(input$run_sc, {
    vec <- c()
    if (!is.null(input$file_sc)) {
      d <- tryCatch({
        read.table(input$file_sc$datapath, header = FALSE, sep = "", 
                   stringsAsFactors = FALSE, fill = TRUE, quote = "\"'")
      }, error = function(e) {
        tryCatch({ read.csv(input$file_sc$datapath, header = FALSE, stringsAsFactors = FALSE) }, 
                 error = function(e2) NULL)
      })
      
      if (!is.null(d) && nrow(d) > 0) {
        if (ncol(d) == 1) {
          vec <- as.character(d[[1]])
        } else if (ncol(d) >= 2) {
          val2 <- suppressWarnings(as.numeric(d[1, 2]))
          sr <- if (is.na(val2)) 2 else 1
          sub_d <- d[sr:nrow(d), ]
          názvy <- as.character(sub_d[[1]])
          počty <- suppressWarnings(as.numeric(sub_d[[2]]))
          validné <- !is.na(počty) & počty > 0
          if (any(validné)) { vec <- rep(názvy[validné], počty[validné]) } else { vec <- as.character(as.vector(t(d))) }
        }
      }
    } else if (input$manual_sc != "") {
      vec <- trimws(strsplit(input$manual_sc, ",")[[1]])
    }
    
    vec <- vec[!is.na(vec) & vec != ""]
    if(length(vec) > 0) {
      rv$data <- vec
      rv$analyzed <- TRUE
      rv$t1 <- input$t1_sc
      rv$t2 <- input$t2_sc
    }
  })
  
  output$results_sc <- renderUI({
    if (!rv$analyzed || is.null(rv$data)) return(NULL)
    div(class = "result-box",
        fluidRow(
          column(5,
                 h3("Výsledky Analýzy", style="color: #026112;"),
                 htmlOutput("text_sc"), hr(),
                 h4("Predpoveď nových druhov:"),
                 tableOutput("table_sc"),
                 wellPanel(style = "background: #e8f5e9; border-left: 5px solid #026112;",
                           strong("Odporúčanie:"), textOutput("rec_sc")
                 )
          ),
          column(7, 
                 plotOutput("plot_sc", height = "700px"),
                 # --- DOWNLOAD BOX PRE VEDECKÚ ČASŤ ---
                 div(class = "download-box",
                     h4("Stiahnuť graf", style = "color: #026112;"),
                     downloadButton("dl_sc_graf", "Stiahnuť analýzu (PNG)", 
                                    class = "btn-science", 
                                    style = "width: 100%; color: white;")
                 )
          )
        )
    )
  })
  
  output$text_sc <- renderUI({
    req(rv$data)
    v <- rv$data; aic <- zisti_najlepsi_model_aic(v); params <- vypocitaj_fisher_parametre(v); pomer <- rv$t2 / rv$t1
    fisher_info <- if (!is.na(params$alpha)) {
      if (params$use_k && params$k > 0) {
        paste0("<br><b>Fisherov test:</b> k = ", round(params$k, 4), " (neg. binomický model)", "<br><b>α (Fisher):</b> ", round(params$alpha, 2))
      } else {
        paste0("<br><b>Fisherov test:</b> k ≈ 0 (log-sériový model)", "<br><b>α (Fisher):</b> ", round(params$alpha, 2))
      }
    } else ""
    gt_info <- if (pomer > 1) "<br><b>Good-Toulmin:</b> t > 1, smoothed verzia" else "<br><b>Good-Toulmin:</b> t ≤ 1, klasická verzia"
    HTML(paste0("<b>Počet jedincov (n):</b> ", length(v), "<br><b>Počet druhov (S<sub>obs</sub>):</b> ", length(unique(v)), "<br><b>Pomer extrapolácie:</b> ", round(pomer, 2), "<br><b>Detegovaný model (AIC):</b> ", aic, fisher_info, gt_info))
  })
  
  output$table_sc <- renderTable({
    req(rv$data)
    v <- rv$data
    data.frame(Metoda = c("Fisher", "Good-Toulmin", "Subsampling (Simulácia)", "Subsampling (Analytický)"),
               Odhad = round(c(predikuj_fisher(v, rv$t1, rv$t2), predikuj_gt(v, rv$t1, rv$t2), predikuj_subsampling_weighted(v, rv$t1, rv$t2), predikuj_subsampling_analyticky(v, rv$t1, rv$t2)), 1))
  }, bordered = TRUE, hover = TRUE)
  
  output$rec_sc <- renderText({
    req(rv$data)
    v <- rv$data; pomer <- rv$t2 / rv$t1; aic <- zisti_najlepsi_model_aic(v)
    msg <- paste0("Dáta: ", aic, ". ")
    if(pomer <= 1) msg <- paste(msg, "Pomer t ≤ 1: optimálny Good-Toulmin.") else if(pomer <= 2) msg <- paste(msg, "Pomer 1 < t ≤ 2: Smoothed GT/Fisher.") else msg <- paste(msg, "Pomer t > 2: Fisher najspoľahlivejší.")
    return(msg)
  })
  
  output$plot_sc <- renderPlot({
    req(rv$data)
    vykresli_graf_sc(rv$data)
  })
  
  # --- DOWNLOAD HANDLER PRE VEDECKÚ ČASŤ ---
  output$dl_sc_graf <- downloadHandler(
    filename = function() "vedecka_analyza.png",
    content = function(file) {
      png(file, width = 1400, height = 1000, res = 120)
      vykresli_graf_sc(rv$data)
      dev.off()
    }
  )
  
  # --- SIMULÁCIA (PRIDANÁ NOTIFIKÁCIA A ITERÁCIE) ---
  data_sim <- eventReactive(input$run_sim, {
    S <- input$S_sim; t1 <- input$t1_sim; t2 <- input$t2_sim; iter <- input$iter_sim
    
    # Načítanie parametrov z rozšírených nastavení (s fallback na pôvodné hodnoty)
    p_ls_theta <- if(!is.null(input$param_ls_theta)) input$param_ls_theta else 0.95
    p_ln_meanlog <- if(!is.null(input$param_ln_meanlog)) input$param_ln_meanlog else 3.0
    p_ln_sdlog <- if(!is.null(input$param_ln_sdlog)) input$param_ln_sdlog else 1.0
    p_gm_p <- if(!is.null(input$param_gm_p)) input$param_gm_p else 0.1
    p_bs_N <- if(!is.null(input$param_bs_N)) input$param_bs_N else 40000
    p_yl_rho <- if(!is.null(input$param_yl_rho)) input$param_yl_rho else 1.5
    
    # Notifikácia pre používateľa počas počítania iterácií
    id <- showNotification(paste("Prebiehajú simulácie (", iter, "iterácií pre každý model), čakajte prosím..."), 
                           duration = NULL, type = "message")
    on.exit(removeNotification(id))
    
    s_ls <- generuj_logseriove(S, p_ls_theta)
    s_ln <- generuj_lognormalne(S, p_ln_meanlog, p_ln_sdlog)
    s_gm <- generuj_geometricke(S, p_gm_p)
    s_bs <- generuj_macarthur_stick(S, p_bs_N)
    s_yl <- generuj_yule(S, p_yl_rho)
    
    res <- rbind(
      analyzuj_jadro_simulacie(s_ls, t1, t2, "Log-Séria", iter), 
      analyzuj_jadro_simulacie(s_ln, t1, t2, "Log-Norm", iter), 
      analyzuj_jadro_simulacie(s_gm, t1, t2, "Geom", iter), 
      analyzuj_jadro_simulacie(s_bs, t1, t2, "MacArthur", iter), 
      analyzuj_jadro_simulacie(s_yl, t1, t2, "Yule", iter)
    )
    list(res=res, w=list(s_ls, s_ln, s_gm, s_bs, s_yl), params=list(t1=t1, t2=t2))
  })
  
  output$results_sim <- renderUI({
    req(data_sim())
    div(class = "result-box", 
        fluidRow(
          column(5, 
                 h3("Priemerná Presnosť (po iteráciách)", style="color: #2980b9;"), 
                 div(class = "table-container", tableOutput("table_sim")), 
                 plotOutput("plot_err_sim", height="300px")
          ), 
          column(7, 
                 h3("Vygenerované Svety (Ukážka 1 behu)", style="color: #2980b9;"), 
                 plotOutput("plot_worlds_sim", height="800px"),
                 # --- DOWNLOAD BOX PRE SIMULAČNÚ ČASŤ ---
                 div(class = "download-box",
                     h4("Stiahnuť grafy", style = "color: #2980b9;"),
                     fluidRow(
                       column(4, downloadButton("dl_logseries", "Log-Séria", class = "btn-sim", style = "width: 100%; margin-bottom: 8px; color: white;")),
                       column(4, downloadButton("dl_lognorm", "Log-Norm", class = "btn-sim", style = "width: 100%; margin-bottom: 8px; color: white;")),
                       column(4, downloadButton("dl_geom", "Geometrické", class = "btn-sim", style = "width: 100%; margin-bottom: 8px; color: white;"))
                     ),
                     fluidRow(
                       column(4, downloadButton("dl_macarthur", "MacArthur", class = "btn-sim", style = "width: 100%; margin-bottom: 8px; color: white;")),
                       column(4, downloadButton("dl_yule", "Yule", class = "btn-sim", style = "width: 100%; margin-bottom: 8px; color: white;")),
                       column(4, downloadButton("dl_all", "Všetky (ZIP)", class = "btn-darkblue", style = "width: 100%; margin-bottom: 8px; color: white;"))
                     )
                 )
          )
        )
    )
  })
  
  output$table_sim <- renderTable({ 
    req(data_sim()) 
    df <- data_sim()$res[, c("Model", "Pred_Fisher", "Pred_GT", "Pred_Sub", "Pred_Sub_Analyticky", "REALITA")]
    colnames(df) <- c("Model", "Priem. Fisher", "Priem. GT", "Priem. Sub", "Priem. Sub(A)", "Priem. Realita")
    df 
  }, bordered = TRUE, spacing = 'xs')
  
  output$plot_err_sim <- renderPlot({ 
    req(data_sim())
    res <- data_sim()$res
    cols <- c("#0dc4fc", "#FF8C00", "#2802e3", "#8A2BE2")
    errs <- t(res[,c("Chyba_Fisher", "Chyba_GT", "Chyba_Sub", "Chyba_Sub_Analyticky")])
    colnames(errs) <- res$Model
    errs[abs(errs)>100] <- 100*sign(errs[abs(errs)>100])
    
    y_lims <- c(min(0, min(errs, na.rm=TRUE)) - 20, max(0, max(errs, na.rm=TRUE)) + 20)
    
    bp <- barplot(errs, beside=TRUE, col=cols, border=NA, main="Priemerné chyby odhadu (%)", las=2, ylim=y_lims)
    text(x = bp, y = errs, labels = round(errs, 1), pos = ifelse(errs >= 0, 3, 1), cex = 0.7, offset = 0.3, xpd = TRUE)
    
    legend("topright", legend=c("Fish","GT","Sub","Sub(A)"), fill=cols, cex=0.7, bty="n")
    abline(h=0) 
  })
  
  output$plot_worlds_sim <- renderPlot({
    req(data_sim()); w <- data_sim()$w; pars <- data_sim()$params; t1 <- pars$t1; t2 <- pars$t2; rate <- 0.1
    par(mfrow=c(5,2), mar=c(4,4,2,1))
    plot_w <- function(p, t) {
      cnt <- sort(table(p), decreasing=T); sad <- table(cnt); x<-as.numeric(names(sad)); y<-as.numeric(sad); df<-data.frame(x,y); df<-df[order(df$x),]; df<-df[df$x<=20,]; if(nrow(df)==0) df<-data.frame(x=1,y=0)
      barplot(df$y, names.arg=df$x, col="#3498db", border=NA, 
              main=paste("Distribúcia abundancií modelu", t),
              xlab="Počet jedincov (n)", ylab="Počet druhov (S_n)")
      N_tot <- length(p); n1 <- round(N_tot * rate * t1); n2 <- round(N_tot * rate * t2); if(n1+n2 > N_tot) { scale<-N_tot/(n1+n2); n1<-floor(n1*scale); n2<-floor(n2*scale) }
      cely_lov <- sample(p); je_novy <- !duplicated(cely_lov); kumulativny_pocet <- cumsum(je_novy)
      plot(1:length(kumulativny_pocet), kumulativny_pocet, type="l", col="gray", lwd=1, 
           main=paste("Krivka objavovania modelu", t), 
           xlab="Počet jedincov", ylab="Kumulatívny počet druhov")
      if (n1 > 0) lines(1:n1, kumulativny_pocet[1:n1], col="#0dc4fc", lwd=2)
      if (n2 > 0 && (n1+n2) <= length(kumulativny_pocet)) lines((n1+1):(n1+n2), kumulativny_pocet[(n1+1):(n1+n2)], col="#e67e22", lwd=2)
      abline(v=n1, lty=2, col="black"); legend("bottomright", legend=c("Tréning", "Test", "Zvyšok"), fill=c("#0dc4fc","#e67e22","gray"), cex=0.6, bty="n")
    }
    plot_w(w[[1]], "Log-Séria"); plot_w(w[[2]], "Log-Norm"); plot_w(w[[3]], "Geometrické"); plot_w(w[[4]], "MacArthur"); plot_w(w[[5]], "Yule")
  })
  
  # --- DOWNLOAD HANDLERS PRE SIMULAČNÚ ČASŤ ---
  uloz_graf_sim <- function(file, populacia, nazov, t1, t2) {
    png(file, width = 1600, height = 500, res = 120)
    vykresli_model_sim(populacia, nazov, t1, t2)
    dev.off()
  }
  
  output$dl_logseries <- downloadHandler(
    filename = function() "logseries_combined.png",
    content = function(file) {
      d <- data_sim(); uloz_graf_sim(file, d$w[[1]], "Log-Séria", d$params$t1, d$params$t2)
    }
  )
  output$dl_lognorm <- downloadHandler(
    filename = function() "lognorm_combined.png",
    content = function(file) {
      d <- data_sim(); uloz_graf_sim(file, d$w[[2]], "Log-Norm", d$params$t1, d$params$t2)
    }
  )
  output$dl_geom <- downloadHandler(
    filename = function() "geom_combined.png",
    content = function(file) {
      d <- data_sim(); uloz_graf_sim(file, d$w[[3]], "Geometrické", d$params$t1, d$params$t2)
    }
  )
  output$dl_macarthur <- downloadHandler(
    filename = function() "brokenstick_combined.png",
    content = function(file) {
      d <- data_sim(); uloz_graf_sim(file, d$w[[4]], "MacArthur", d$params$t1, d$params$t2)
    }
  )
  output$dl_yule <- downloadHandler(
    filename = function() "yule_combined.png",
    content = function(file) {
      d <- data_sim(); uloz_graf_sim(file, d$w[[5]], "Yule", d$params$t1, d$params$t2)
    }
  )
  
  # --- ZIP so všetkými 5 grafmi naraz ---
  output$dl_all <- downloadHandler(
    filename = function() "vsetky_grafy.zip",
    content = function(file) {
      d <- data_sim()
      
      # Vytvoríme čerstvý dočasný priečinok len pre tento ZIP
      tmpdir <- file.path(tempdir(), paste0("zip_", as.integer(Sys.time())))
      dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
      on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
      
      nazvy <- c("logseries_combined.png", "lognorm_combined.png", "geom_combined.png",
                 "brokenstick_combined.png", "yule_combined.png")
      modely <- c("Log-Séria", "Log-Norm", "Geometrické", "MacArthur", "Yule")
      cesty <- file.path(tmpdir, nazvy)
      
      # Vygenerujeme všetkých 5 PNG súborov
      for (i in 1:5) {
        uloz_graf_sim(cesty[i], d$w[[i]], modely[i], d$params$t1, d$params$t2)
      }
      
      # Použijeme balíček `zip` (čistá R implementácia, funguje všade)
      if (requireNamespace("zip", quietly = TRUE)) {
        zip::zipr(zipfile = file, files = cesty)
      } else {
        # Fallback na base zip(), ak by používateľ nemal balíček `zip`
        old_wd <- setwd(tmpdir)
        on.exit(setwd(old_wd), add = TRUE)
        utils::zip(zipfile = file, files = nazvy)
      }
    }
  )
}

shinyApp(ui, server)