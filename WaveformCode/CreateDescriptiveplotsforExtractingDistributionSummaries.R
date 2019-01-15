p1 <- BC_PlotCreateRRTimesPlots( RPeaksStruct = RPeakData, MetaData = MetaData)
p1 <- p1 + geom_vline(xintercept = RRDistributionSummariesOutput[13 , 1] , color = 'blue') +
  geom_vline(xintercept = RRDistributionSummariesOutput[13 , 25] , color = 'red')

p2 <- ggplot(data.frame( RRTimes =tmpKde$eval.points , denisty = tmpKde$estimate) , aes(RRTimes , denisty)) + geom_line( color = 'blue') + ggtitle('Far from AF')
p3 <- ggplot(data.frame( RRTimes =tmpKde$eval.points , denisty = tmpKde$estimate) , aes(RRTimes , denisty)) + geom_line( color = 'red') + ggtitle('Close to AF')
p4 <- ggplot(data.frame( RRTimes =tmpKde$eval.points , denisty = tmpKde$estimate) , aes(RRTimes , denisty))+ geom_line( color = 'purple') + ggtitle('In AF')


x11(20,14)
lay <- rbind(c(1,1,1),
             c(2,3,4))
grid.arrange(p1 , p2+ xlim(0.25,1.5) + ylim(0,12.5), p3+ xlim(0.25,1.5)+ ylim(0,12.5),p4+ xlim(0.25,1.5)+ ylim(0,12.5), layout_matrix = lay)
