# Logistic-regression
Program to perform Logistic regression
Logistic regression 
 
    INPUT: n (number of data point, D), mx1, vx1, my1, vy1, mx2, vx2, my2, vy2 (m: mean, v: variance)
 
    FUNCTION: 
                   a. Generate n data point: D1= {(x1, y1), (x2 ,y2), ..., (xn, yn) }, where x and y are independently sampled from N(mx1, vx1) and N(my1, vy1) respectively. (use the Gaussian random number generator you did for homework 3.).
                   b. Generate n data point: D2= {(x1, y1), (x2 ,y2), ..., (xn, yn) }, where x and y are independently sampled from N(mx2, vx2) and N(my2, vy2) respectively. 
                   c. Use Logistic regression to separate D1 and D2. You should implement both Newton's and steepest gradient descent method during optimization, i.e., when the Hessian is singular, use steepest descent for instead. You should come up with a reasonable rule to determine convergence. (a simple run out of the loop should be used as the ultimatum) 
 
     OUTPUT: The confusion matrix and the sensitivity and specificity of the logistic regression applied to the training data D.
