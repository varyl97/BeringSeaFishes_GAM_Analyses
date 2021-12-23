euclidean.distance<-function(a,b) sqrt(sum((a-b)^2))


euclidean.distance<-function(start.sal,start.temp,
                             end.sal,end.temp){
  delta.sal<- start.sal-end.sal
  delta.temp<- start.temp-end.temp
  distance<- sqrt(sum(delta.temp-delta.sal)^2)
  distance
}