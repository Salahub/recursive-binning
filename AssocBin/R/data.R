##' De-Garched S&P 500 returns
##'
##' This data is the result of code from the 'zenplots' package to
##' process S&P 500 consituent stock returns into uniform
##' pseudo-observations for measuring association.
##'
##' @format
##' A matrix with 755 rows and 461 columns, the rows correspond to
##' dates between 2007 and 2009 and the columns correspond to the
##' different S&P 500 constituent stocks.
##'
##' @usage data(sp500pseudo)
'sp500pseudo'

##' Heart Disease Diagnosis Data
##'
##' This data (adapted from the UCI Machine Learning Repository at
##' https://archive.ics.uci.edu/) presents a single data frame
##' reporting heart disease diagnosis results for patients from
##' studies carried out by Andras Janosi at the Hungarian
##' Institute of Cardiology; William Steinbrunn and Matthias
##' Pfisterer at the University Hospitals of Zurich and Basel;
##' and two separate studies by Robert Detrano carried out at the
##' Cleveland Clinic Foundation and Long Beach V.A. Medical Center.
##' The data contains measurements of 15 variables collected on 920
##' participants:
##' \describe{
##'    \item{age}{Age in years}
##'    \item{sex}{Sex}
##'    \item{cp}{Reported chest pain type: typical angina, non-typical angina, non-angina, or no pain}
##'    \item{trestbps}{Resting blood pressure (mmHg on admission to hospital)}
##'    \item{chol}{Serum cholesterol in mg/dl}
##'    \item{fbs}{Indicator of fasting blood sugar >120 mg/dl}
##'    \item{restecg}{Resting electrocardiographic results: normal, indicating ventricular hypertrophy, or displaying ST-T wave abnormality}
##'    \item{thalach}{Maximum measured heart rate}
##'    \item{exang}{Indicator of exercise induced angina}
##'    \item{oldpeak}{ST wave depression induced by exercise relative to rest}
##'    \item{slope}{The slope of the ST segment during peak exercise}
##'    \item{ca}{Number of major blood vessels coloured by fluoroscopy}
##'    \item{thal}{Type of heart defect}
##'    \item{num}{Diagnosis of heart disease. Values greater than one indicate heart disease of different sorts while a value of zero indicates no heart disease}
##'    \item{study}{The study where the participant's data was collected}
##' }
##'
##' @format
##' A matrix with 920 rows and 15 columns, with each row reporting
##' measurements for a participant in one of the heart disease
##' studies.
##'
##' @usage data(heart)
'heart'
