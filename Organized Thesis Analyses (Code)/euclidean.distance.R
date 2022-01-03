euclidean.distance<-function(start.sal,start.temp,
                             end.sal,end.temp){
  delta.sal<- start.sal-end.sal
  delta.temp<- start.temp-end.temp
  distance<- sqrt(delta.temp^2+delta.sal^2)
  distance
}