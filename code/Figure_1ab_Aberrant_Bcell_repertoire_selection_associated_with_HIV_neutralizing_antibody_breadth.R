
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(argparser, quietly=TRUE)

p <- arg_parser("Pre process data and create figures highlighting mutation levels for subjects in HIV study cohort")
p <- add_argument(p, "subject_mutationleve_data", help="csv file containing per subject V segment mutation levels")
p <- add_argument(p, "meta_data", help="csv file containing information on HIV status of subjects in cohort")

argv <- parse_args(p)

# read in required files i.e mutation level data of subjects and metadata containing disease phenotype
mutation_level = fread(data.dir = argv$subject_mutationleve_data, sep =',')
metadata = read.csv(data.dir = argv$meta_data, header=T)

#remove unwanted string from subject ID's in mutation level file
mutation_level$subject = gsub('boydlab:', "", mutation_level$subject)

#merge meta data to mutation level data
mutation_level_meta = merge(mutation_level, metadata, by ='subject')

#make groups based on subjects that are HIV negative, noNab and BNabs
HIV_neagtive = subset(mutation_level_meta, diagnosis =='HIV Negative')
noNab = subset(mutation_level_meta, diagnosis =='HIV Non Neutralizing')
BNab = subset(mutation_level_meta, diagnosis =='HIV Broad Neutralizing')

HIV_neagtive = HIV_neagtive %>% mutate(diagnosis2 = 'HIV_neg' )
noNab =noNab %>% mutate(diagnosis2 = 'noNab' )   
BNab= BNab %>% mutate(diagnosis2 = 'BNab' ) 

subject_mutation_data= rbind(HIV_neagtive, noNab, BNab)

#Remove unwanted characters from names of  isotypes
subject_mutation_data$isotype = ifelse(subject_mutation_data$isotype =='IGHG1', 'IgG1',
  ifelse(subject_mutation_data$isotype =='IGHG2', 'IgG2',
    ifelse(subject_mutation_data$isotype =='IGHG3', 'IgG3',
      ifelse(subject_mutation_data$isotype =='IGHG4', 'IgG4',
        ifelse(subject_mutation_data$isotype =='IGHD', 'IgD',
          ifelse(subject_mutation_data$isotype =='IGHM', 'IgM',
            ifelse(subject_mutation_data$isotype =='IGHE', 'IgE',
              ifelse(subject_mutation_data$isotype =='IGHA1', 'IgA1',
                ifelse(subject_mutation_data$isotype =='IGHA2', 'IgA2',  
                  subject_mutation_data$isotype)))))))))


#change the levels for subject mutation level diagnosis2 column
subject_mutation_data$diagnosis2 = factor(subject_mutation_data$diagnosis2, levels = c('HIV_neg', 'noNab', 'BNab'))

#subset subject_mutation_data based on isotypes
#Naive isotypes
subject_mutation_data_naive = subset(subject_mutation_data, isotype == 'IgM'|isotype == 'IgD')
#IgA isotype
subject_mutation_data_IgA   =subset(subject_mutation_data, isotype == 'IgA1'|isotype == 'IgA2')
#IgG isotype
subject_mutation_data_IgG   =subset(subject_mutation_data, isotype == 'IgG1'|isotype == 'IgG2'|isotype == 'IgG3')
#IgE isotype
subject_mutation_data_IgE   =subset(subject_mutation_data, isotype == 'IgE')


#create boxplots, figure 1a
col = c( noNab= "#009E73", HIV_neg="#D55E00", BNab="#56B4E9")
naive = ggplot(subject_mutation_data_naive, aes(x=diagnosis2, y=mutation_level.median.mean, fill = diagnosis2))+
  geom_boxplot(outlier.shape = NA, width =1)+ 
  geom_point( position = position_jitter( width =0.09,height =0 ), size =0.3)+
  theme_bw(base_size=8)+
  theme(axis.text=element_text(size=7, face= 'bold'),axis.title=element_text(size=7),
        legend.text = element_text(size=7, face ='bold'), legend.title = element_text(size=7, face = 'bold'))+
  ylab('Average SHM frequency (%)')  +
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')+
  facet_grid(. ~ isotype, scales = "free", space='free') +
  scale_fill_manual(values= col)

IgA=ggplot(subject_mutation_data_IgA, aes(x=diagnosis2, y=mutation_level.median.mean, fill = diagnosis2))+
  geom_boxplot(outlier.shape = NA, width =1)+ 
  geom_point( position = position_jitter( width =0.09,height =0 ), size =0.3)+
  theme_bw(base_size=8)+
  theme(axis.text=element_text(size=7, face= 'bold'),axis.title=element_text(size=7),
        legend.text = element_text(size=7, face ='bold'), legend.title = element_text(size=7, face = 'bold'))+
  ylab('Average SHM frequency (%)')  +
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')+
  facet_grid(. ~ isotype, scales = "free", space='free') +
  scale_fill_manual(values= col)

my_comparisons <- list( c("HIV_neg", "noNab"), c("noNab", "BNab"), c("HIV_neg", "BNab") )
IgG=ggplot(subject_mutation_data_IgG, aes(x=diagnosis2, y=mutation_level.median.mean, fill = diagnosis2))+
  geom_boxplot(outlier.shape = NA, width =1)+ 
  geom_point( position = position_jitter( width =0.09,height =0 ), size =0.3)+
  theme_bw(base_size=8)+
  theme(axis.text=element_text(size=7, face= 'bold'),axis.title=element_text(size=7),
        legend.text = element_text(size=7, face ='bold'), legend.title = element_text(size=7, face = 'bold'))+
  ylab('Average SHM frequency (%)')  +
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = scales::percent, n.breaks=5)+
  theme(legend.position = 'none')+
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)+
  facet_grid(. ~ isotype, scales = "free", space='free') +
  scale_fill_manual(values= col)

IgE=ggplot(subject_mutation_data_IgE, aes(x=diagnosis2, y=mutation_level.median.mean, fill = diagnosis2))+
  geom_boxplot(outlier.shape = NA, width =1)+ 
  geom_point( position = position_jitter( width =0.09,height =0 ), size =0.3)+
  theme_bw(base_size=8)+
  theme(axis.text=element_text(size=7, face= 'bold'),axis.title=element_text(size=7),
        legend.text = element_text(size=7, face ='bold'), legend.title = element_text(size=7, face = 'bold'))+
  ylab('Average SHM frequency (%)')  +
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')+
  facet_grid(. ~ isotype, scales = "free", space='free') +
  scale_fill_manual(values= col)

## plot 99th percentile of somatic hypermutation per subject, figure 1b
my_comparisons <- list( c("HIV_neg", "noNab"), c("noNab", "BNab"), c("HIV_neg", "BNab") )
IgG_b=ggplot(subject_mutation_data_IgG, aes(x=diagnosis2, y=mutation_level.median.percentie99, fill = diagnosis2))+
  geom_boxplot(outlier.shape = NA, width =1)+ 
  geom_point( position = position_jitter( width =0.09,height =0 ), size =0.3)+
  theme_bw(base_size=8)+
  theme(axis.text=element_text(size=7, face= 'bold'),axis.title=element_text(size=7),
        legend.text = element_text(size=7, face ='bold'), legend.title = element_text(size=7, face = 'bold'))+
  ylab('Average SHM frequency (%)')  +
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')+
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)+
  facet_grid(. ~ isotype, scales = "free", space='free') +
  scale_fill_manual(values= col)

naive_b=ggplot(subject_mutation_data_naive, aes(x=diagnosis2, y=mutation_level.median.percentie99, fill = diagnosis2))+
  geom_boxplot(outlier.shape = NA, width =1)+ 
  geom_point( position = position_jitter( width =0.09,height =0 ), size =0.3)+
  theme_bw(base_size=8)+
  theme(axis.text=element_text(size=7, face= 'bold'),axis.title=element_text(size=7),
        legend.text = element_text(size=7, face ='bold'), legend.title = element_text(size=7, face = 'bold'))+
  ylab('Average SHM frequency (%)')  +
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')+
  facet_grid(. ~ isotype, scales = "free", space='free') +
  scale_fill_manual(values= col)

IgA_b=ggplot(subject_mutation_data_IgA, aes(x=diagnosis2, y=mutation_level.median.percentie99, fill = diagnosis2))+
  geom_boxplot(outlier.shape = NA, width =1)+ 
  geom_point( position = position_jitter( width =0.09,height =0 ), size =0.3)+
  theme_bw(base_size=8)+
  theme(axis.text=element_text(size=7, face= 'bold'),axis.title=element_text(size=7),
        legend.text = element_text(size=7, face ='bold'), legend.title = element_text(size=7, face = 'bold'))+
  ylab('Average SHM frequency (%)')  +
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')+
  facet_grid(. ~ isotype, scales = "free", space='free') +
  scale_fill_manual(values= col)

IgE_b=ggplot(subject_mutation_data_IgE, aes(x=diagnosis2, y=mutation_level.median.percentie99, fill = diagnosis2))+
  geom_boxplot(outlier.shape = NA, width =1)+ 
  geom_point( position = position_jitter( width =0.09,height =0 ), size =0.3)+
  theme_bw(base_size=8)+
  theme(axis.text=element_text(size=7, face= 'bold'),axis.title=element_text(size=7),
        legend.text = element_text(size=7, face ='bold'), legend.title = element_text(size=7, face = 'bold'))+
  ylab('Average SHM frequency (%)')  +
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')+
  facet_grid(. ~ isotype, scales = "free", space='free') +
  scale_fill_manual(values= col)

# save figures to local folder
