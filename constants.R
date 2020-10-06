source('~/gnomad_lof/R/constants.R')
library(ggthemes)
# BiocManager::install('ggbio')

theme_set(theme_classic())

output_type = function(output_format='png', fname, height, width, res=300, units='in', onefile=TRUE) {
  if (output_format == 'png') {
    return(png(fname, height=height, width=width, res=res, units=units))
  } else {
    return(pdf(fname, height=height, width=width, onefile=onefile))
  }
}

pops = c('AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID')
ukb_gnomad_pop_mapping = c('AFR' = 'afr', 'AMR' = 'amr', 'CSA' = 'sas', 'EAS' = 'eas', 'EUR' = 'nfe', 'MID' = 'mde')
pops_with_gnomad_data = c('AFR', 'AMR', 'EAS', 'EUR')
ukb_pop_colors = pop_colors[ukb_gnomad_pop_mapping]
names(ukb_pop_colors) = pops
ukb_pop_colors['MID'] = '#EEA9B8'

pop_color_scale = scale_color_manual(values=ukb_pop_colors, name='Population')
pop_fill_scale = scale_fill_manual(values=ukb_pop_colors, name='Population')
xyline = geom_abline(slope = 1, intercept = 0, linetype='dashed')

trait_types = c("biomarkers", "continuous", "categorical", "phecode", "icd10", "prescriptions")
trait_type_colors = few_pal('Dark')(6)
names(trait_type_colors) = trait_types
trait_type_names = c('Biomarkers', 'Continuous', 'Categorical', 'Disease (phecode)', 'Disease (ICD)', 'Prescriptions')
names(trait_type_names) = trait_types
trait_type_colors['hide'] = 'gray'
trait_color_scale = scale_color_manual(breaks = trait_types, values=trait_type_colors, name='Trait type', labels=trait_type_names)
trait_fill_scale = scale_fill_manual(breaks = trait_types, values=trait_type_colors, name='Trait type', labels=trait_type_names)

loglog_breaks = c(0:10, 20, 50, 100, 200, 400)
ll_to_pvalue = function(yt) {
  return(if_else(yt >= 10, 10 ^ (yt / 10), yt))
}
pvalue_to_ll = function(y) {
  return(if_else(y >= 10, 10*log10(y), y))
}
gwas_loglog_trans = function() {
  scales::trans_new("gwas_loglog", transform = pvalue_to_ll, inverse = ll_to_pvalue)
}

get_ukb_data_url = function() {
  return(paste0('https://storage.googleapis.com/ukb-diverse-pops-public-free/sumstats_qc_analysis/'))
}

check_size = function(local_fname, url) {
  if (is.character(getURL("www.google.com"))) {
    response = httr::HEAD(url)
    remote_size = httr::headers(response)[["Content-Length"]]
    local_size = file.size(local_fname)
    same_size = remote_size == local_size
    if (same_size) {
      inform(paste0(local_fname, ' is the same size as remote file (', remote_size, '). Using cached version...'))
    } else {
      inform(paste0(local_fname, ' is a different size (', local_size, ') from remote file (', remote_size, '). Downloading new version...'))
    }
    return(same_size)
  } else {
    inform(paste('No internet connection. Using cached version of', local_fname))
  }
}

get_or_download_ukb_file = function(base_fname, subfolder='', local_name='', use_local=F) {
  fname = paste0(data_dir, ifelse(local_name != '', local_name, base_fname))
  url = paste0(get_ukb_data_url(), subfolder, base_fname)
  if (!file.exists(fname)) {
    download.file(url, fname)
  } else {
    if (!use_local) {
      if (!check_size(fname, url)) {
        download.file(url, fname)
      }
    }
  }
  return(fname)
}

load_ukb_file = function(base_fname, subfolder='', local_name='', force_cols=cols(), use_local=F) {
  fname = get_or_download_ukb_file(base_fname, subfolder, local_name, use_local=use_local)
  if (endsWith(fname, '.gz') | endsWith(fname, '.bgz')) {
    fname = gzfile(fname)
  }
  data = read_delim(fname, delim='\t', col_types = force_cols)
  return(data)
}

blue = '#002F6C'
low_color = '#EEEEEE'
high_color = blue

pretty_axis_format = function(x) {
  return(case_when(!is.finite(x) | (x < 1e-4) ~ scientific(x),
                   x >= 1e-1 ~ percent(x, accuracy = 1e+1),
                   x >= 1e-2 ~ percent(x, accuracy = 1e+0),
                   x >= 1e-3 ~ percent(x, accuracy = 1e-1),
                   x >= 1e-4 ~ percent(x, accuracy = 1e-2),
                   TRUE ~ 'blah' # percent(x, accuracy = 1e+2)
  ))
}

pretty_axis_format_binned = function(x) {
  x = x * (1 + 9 * corrected)
  return(case_when(
    x == 1e-4 ~ '0.001-0.01%',
    x == 1e-3 ~ '0.01-0.1%',
    x == 1e-2 ~ '0.1-1%',
    x == 1e-1 ~ '1-10%',
    x == 1e-0 ~ '10-100%',
    TRUE ~ 'blah'
  ))
}
pretty_axis_format_binned_corrected = function(x) {
  x = x * 10
  return(case_when(
    x == 1e-4 ~ '0.001-0.01%',
    x == 1e-3 ~ '0.01-0.1%',
    x == 1e-2 ~ '0.1-1%',
    x == 1e-1 ~ '1-10%',
    x == 1e-0 ~ '10-100%',
    TRUE ~ 'blah'
  ))
}

get_color = function(ukb_pop) {
  if (ukb_pop %in% names(ukb_gnomad_pop_mapping)) {
    return(pop_colors[ukb_gnomad_pop_mapping[[ukb_pop]]])
  } else {
    return('gray')
  }
}