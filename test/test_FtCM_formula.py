
from mut_framework.BTagSF_FtCM.FtCM_formula import multiset, submultiset
from mut_framework.BTagSF_FtCM.FtCM_formula import prob_component, prob_submultiset

print "multiset(1,2) = " + str(multiset(1,2))
print "multiset(2,3) = " + str(multiset(2,3))
print "multiset(3,3) = " + str(multiset(3,3))
print "multiset(3,1) = " + str(multiset(3,1))
print "submultiset((1,1,1),3) = " + str(submultiset((1,1,1),3))
print "submultiset((3,3,3),3) = " + str(submultiset((3,3,3),3))
print "submultiset((1,1,2),1) = " + str(submultiset((1,1,2),1))

print "prob_component((1,0,0),(1,1,2)) = " +str(prob_component((1,0,0),(1,1,2))) 
print "prob_component((0,1,0),(1,1,2)) = " +str(prob_component((0,1,0),(1,1,2))) 
print "prob_component((0,0,1),(1,1,2)) = " +str(prob_component((0,0,1),(1,1,2))) 
print "prob_submultiset((1,1,2), 1) = " +str(prob_submultiset((1,1,2), 1)) 
print "prob_submultiset((1,1,2), 1).expand() = " + str(prob_submultiset((1,1,2), 1).expand()) 
