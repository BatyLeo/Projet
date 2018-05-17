def RIEMANN(lambda1,state1,lambda2,state2):
    lambda0=(lambda1+lambda2)/2. ## coefficient = wave-speed
    rightf=max(lambda0,0)*state1+min(lambda0,0)*state2 - lambda0*state2
    leftf=-(max(lambda0,0)*state1+min(lambda0,0)*state2 - lambda0*state1 )
    return [leftf,rightf,lambda0] ## end of Lax-Friedrichs Riemann solver
    