
######################################################################
#Check distance of two points in a rectangle of w*l
#The distance should be periodically defined: meaning
#if the point near the opposite edges of the window 
#should also be checked in terms of center distance
#By using different Dist function, we can obmit the edge effect. The distribution would be different.
Dist_no_cyclic<-function(p1,p2,h,w)
{
  return(dist0<-sqrt(sum((p1-p2)^2)))
}