sta_lat <- function(m, se, lab, met, cap = "", minmax = 'min', sci = FALSE, dig = 3, filNam){
  # Dimensions
  nMet <- dim(m)[1]
  nSta <- dim(m)[2]
  sta <- c()
  sta[1] <- "\\begin{table} {\\small"
  sta[2] <- paste0("\\caption{\\enspace ", cap, "}")
  sta[3] <- "    \\begin{center}"
  numCol <- paste0(rep("c", nSta), collapse = "")
  sta[4] <- paste0("        \\begin{tabular*}{\\hsize}{@{\\extracolsep{\\fill}}|l|", numCol, "|}")
  sta[5] <- "            \\hline"
  sta[6] <- "            \\\\[-7pt]"
  sta[7] <- "            \\multicolumn{1}{c}{\\it Method} &"
  for(i in 1:(nSta - 1)){
    sta[7 + i] <- paste0("            \\multicolumn{1}{c}{\\it ", lab[i],"} &")
  }
  sta[7 + nSta] <- paste0("            \\multicolumn{1}{c}{\\it ", lab[nSta],"} \\\\")
  sta[7 + nSta + 1] <- "            \\hline"
  sta[7 + nSta + 2] <- "            \\\\[-5pt]"
  # Obtains the best performer
  bes <- numeric(nSta)
  for(j in 1:nSta){
    if(minmax[j] == 'max'){
      bes[j] <- which.max(m[ ,j])
    } else {
      bes[j] <- which.min(m[ ,j])
    }
  }
  for(i in 1:nMet){
    lin <- paste0(met[i], " &") 
    for(j in 1:(nSta - 1)){
      if(bes[j] != i){
        lin <- paste0(lin, " $", format(m[i, j], scientific = sci, digits = dig), "$ &")
      } else {
        lin <- paste0(lin, " $\\mathbf{", format(m[i, j], scientific = sci, digits = dig), "}$ &")
      }
    }
    if(bes[nSta] != i){
      lin <- paste0(lin, " $", format(m[i, nSta], scientific = sci, digits = dig), "$ \\\\[-5pt]")
    } else {
      lin <- paste0(lin, " $\\mathbf{", format(m[i, nSta], scientific = sci, digits = dig), "}$ \\\\[-5pt]")
    }
    sta[7 + nSta + 2 + i * 2 - 1] <- lin
    lin <- paste0(" ", " &") 
    for(j in 1:(nSta - 1)){
      lin <- paste0(lin, " $_{(", format(se[i, j], scientific = sci, digits = dig), ")}$ &")
    }
    lin <- paste0(lin, " $_{(", format(se[i, nSta], scientific = sci, digits = dig), ")}$ \\\\")
    sta[7 + nSta + 2 + i * 2] <- lin
  }
  sta[7 + nSta + 2 + 2 * nMet + 1] <- "            \\hline"
  sta[7 + nSta + 2 + 2 * nMet + 2] <- "        \\end{tabular*}"
  sta[7 + nSta + 2 + 2 * nMet + 3] <- "    \\end{center}"
  sta[7 + nSta + 2 + 2 * nMet + 4] <- "}\\end{table}"
  
  # Saves to file
  write.table(sta, file = filNam, quote = FALSE, row.names = FALSE)
}