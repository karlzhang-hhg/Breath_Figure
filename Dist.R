
######################################################################
#Check distance of two points in a rectangle of w*l
#The distance should be periodically defined: meaning
#if the point near the opposite edges of the window 
#should also be checked in terms of center distance
#By using different Dist function, we can obmit the edge effect. The distribution would be different.
Dist<-function(p1,p2,h,w)
{
  dist0<-sqrt(sum((p1-p2)^2))
  distw<-sqrt(sum((c(w,0)-abs(p1-p2))^2))
  disth<-sqrt(sum((c(0,h)-abs(p1-p2))^2))
  if (dist0<distw & dist0<disth) return(c(0,0,dist0))
  else
  {
    if(distw<disth)
    {
      return(c(1,0,distw))
    }
    else
    {
      return(c(0,1,disth))
    }
  }
}