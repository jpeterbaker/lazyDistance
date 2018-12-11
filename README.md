Calculates an approximate distance matrix while measuring as few actual distances as possible.

Here's the scenario:

    The distance d between each pair of objects (points or nodes) must be estimated
    
    Each distance can be obtained exactly, but it's very expensive
    
    The "distance" is a proper metric

This code uses metric properties to minimize the number of exact measurements

while guaranteeing that every distance is known within some tolerance

I think method needs two inputs besides the points and distance function:

    1) tolerance
    
    2) how the tolerance is to be applied
    
Each pair of points will have an upper and lower distance bound

    say lo <= d(A,B) <= hi
    
The obvious way to use the tolerance are to iterate until

    criteria(A,B) := (hi-lo)/Z < tol
    
where Z could be chosen different ways

    edgewise: Z = hi
    
    globally: Z = longest measured known edge
    
    safe globally: Z = lower bound on longest actual edge (largest lo value)

If the distance from A-B has been calculated, criteria(A,B) = 0

I suggest that whichever pair has the largest criteria

should have their distance calculated next

Each pair should have an upper and a lower bound

These go out of date as new edges are measured

Here's how the bounds are updated

Say we know ABlo <= AB = <= ABhi (notation generalizes)

Say that CD has just been measured or its bounds have been modified

1) ABlo update

    If CD is long ("forcing" AB apart)
    
        CD <= CA + AB + BD
        
        AB >= CD - AC - BD
        
           >= CDlo - AChi - BDhi
           
        So ABlo can be updated if this quantity is larger than current value
        
        (likewise with C and D reversed)
        
    If CD is short (possibly helping some other long edge to push A,B apart)
    
        CAlo <= CA
        
             <= CD + DB + BA
             
             <= CD + DBhi + BA
             
        AB >= CAlo - CD - DBhi
        
           >= CAlo - CDhi - DBhi
           
        So ABlo can be updated if this quantity is larger than current value
        
        (likewise with C and D reversed)

2) ABhi update

    AB <= AC + CD + BD
    
       <= AChi + CDhi + BDhi
       
    So ABhi can be updated if this quantity is smaller than current value
    
    (likewise with C and D reversed)

