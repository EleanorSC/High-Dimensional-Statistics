# High-Dimensional-Statistics
Intro to high dimensional statistics for the Intermediate Statistics and Data Management 2024 course

[All code for exercises](https://eleanorsc.github.io/High-Dimensional-Statistics/)
## Setup

```r
# Setup your working directory, e.g.
setwd("/Users/username/repos/High_Dimensional_Stats")
```

You will need to have certain dependencies set up for this course. This should not take long to run, but depends on the speed of your computer, your internet connection, and any packages you have installed already - getting this set up before the course is therefore optimal. You’ll need to install R 4.0 or later.

R usually enables package downloads using pre-built binaries. Some times, this is not possible, particularly on Linux and Mac systems. In this case, R package installation often requires additional system dependencies. If you are a Linux user, to ensure that you can download packages using the code below, first run the terminal commands for your distribution here. Note that you will need to use root access (sudo) to install the system dependencies. Mac users may need to use homebrew to install system dependencies, and Windows users may need to install RTools. Ideally, installing packages will proceed without error and you can ignore these steps, but this isn’t always the case.

```r
# Install BiocManager
install.packages("BiocManager")

# Download and read dependencies file
download.file(
  "https://raw.githubusercontent.com/EleanorSC/High-Dimensional-Statistics/main/dependencies.csv",
  destfile = 'dependencies.csv'
)
table <- read.table('dependencies.csv')

# Install dependencies using BiocManager
BiocManager::install(table[[1]])

# Create a directory for data files
dir.create("data", showWarnings = FALSE)

# List of data files to download
data_files <- c(
  "prostate.rds",
  "coefHorvath.rds",
  "methylation.rds"
)

# Download data files into the "data" directory
for (file in data_files) {
  download.file(
    url = file.path(
      "https://raw.githubusercontent.com/EleanorSC/High-Dimensional-Statistics/main/Data",
      file
    ),
    destfile = file.path("data", file)
  )
}
```
- Installs tools (`BiocManager`) for managing biological data analysis packages.
- Downloads a list of necessary R packages from GitHub.
- Installs these R packages.
- Creates a folder (`data`) to store project data.
- Downloads specific data files into this folder from GitHub.
- This code prepares your R environment by setting up all the required packages and data for the project.


## Data Repository

Access the data files used in this project directly from the repository:  [Data](https://github.com/EleanorSC/High-Dimensional-Statistics/tree/main/Data)

This repository contains datasets and associated descriptions used for analyzing high-dimensional statistics with R. Below are the details of the datasets included in this project.

---

## Prostate Cancer Data

[SOURCE](https://github.com/EleanorSC/High-Dimensional-Statistics/tree/main/Data/prostate.rds)
Prostate-specific antigen values and clinical measures for 97 patients hospitalized for a radical prostatectomy. Prostate specimens underwent histological and morphometric analysis.

### **Column Descriptions**
- **lcavol**: log(cancer volume)
- **lweight**: log(prostate weight)
- **age**: age of the patient
- **lbph**: log(benign prostatic hyperplasia amount)
- **svi**: seminal vesicle invasion
- **lcp**: log(capsular penetration)
- **gleason**: Gleason score
- **pgg45**: percentage Gleason scores 4 or 5
- **lpsa**: log(prostate specific antigen)

---

## Methylation Data

[SOURCE](https://github.com/EleanorSC/High-Dimensional-Statistics/tree/main/Data/methylation.rds)
Illumina Human Methylation data from EPIC on sorted peripheral adult blood cell populations. The dataset records DNA methylation assays, which measure the proportion of DNA that carries a methyl mark. The assays are normalized methylation levels (M-values), with:
- Negative values: unmethylated DNA
- Positive values: methylated DNA

### **Data Object**
- **assay(data)**: normalized methylation levels
- **colData(data)**: individual-level metadata

### **Phenotypic Metadata**
- **Sample_Well**: sample well
- **Sample_Name**: name of sample
- **purity**: sample cell purity
- **Sex**: sex of the individual
- **Age**: age in years
- **weight_kg**: weight in kilograms
- **height_m**: height in meters
- **bmi**: body mass index (BMI)
- **bmi_clas**: BMI class
- **Ethnicity_wide**: ethnicity, wide class
- **Ethnic_self**: ethnicity, self-identified
- **smoker**: yes/no indicator of smoker status
- **Array**: type of array from the EPIC array library
- **Slide**: slide identifier

---

## Horvath Data

[SOURCE](https://github.com/EleanorSC/High-Dimensional-Statistics/tree/main/Data/coefHorvath.rds)
Methylation markers across different age groups. This dataset includes CpGmarker variables, which are CpG site encodings.

---

This repository provides resources for advanced statistical modeling, visualization, and interpretation of high-dimensional biological data in R.

---
## License

Distributed under the terms of the [MIT license](LICENSE).
