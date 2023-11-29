#tornado uncertainty intervals in table
library(ggplot2)
library(dplyr)

df<-data.frame(solution_id = c(1,2),
               mean_obj1 = c(2,3),
               lower_obj1 = c(1,2),
               upper_obj1 = c(3,4))%>%
  arrange(mean_obj1)%>%
  mutate(ymin_ref = rank(-mean_obj1)-0.25,
         ymax_ref = rank(-mean_obj1)+0.25)

ggplot()+
  geom_rect(data = df, mapping = aes(xmin = lower_obj1, 
                                     xmax = upper_obj1,
                                     ymin = ymin_ref,
                                     ymax = ymax_ref,
                                     fill = as.factor(solution_id)))+
  geom_point(data = df, aes(x = mean_obj1, y = (ymin_ref+ymax_ref)/2))+
  geom_text(data=df, aes(x = mean_obj1, 
                         y =ymax_ref+.15, 
                         label=paste0(mean_obj1, ' [', lower_obj1, ', ',
                                      upper_obj1, ']'), size=12))+
  geom_hline(yintercept = 1.5)
  
ggplot() + 
  geom_rect(data = tornado_plot_df2,
            aes(ymax=ymax, 
                ymin=ymin, 
                xmax=xmax, 
                xmin=xmin, 
                fill=variable))