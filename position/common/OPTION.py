class POSTIONOPTION():
    '''
    const prcopt_t prcopt_default={ /* defaults processing options */
    PMODE_SINGLE,0,2,SYS_GPS,   /* mode,soltype,nf,navsys */
    15.0*D2R,{{0,0}},           /* elmin,snrmask */
    0,1,1,5,0,10,               /* sateph,modear,glomodear,maxout,minlock,a */
    0,0,0,0,                    /* estion,esttrop,dynamics,tidecorr */
    1,0,0,0,0,                  /* niter,codesmooth,intpref,sbascorr,sbassatsel */
    0,0,                        /* rovpos,refpos */
    {100.0,100.0},              /* eratio[] */
    {100.0,0.003,0.003,0.0,1.0}, /* err[] */
    {30.0,0.03,0.3},            /* std[] */
    {1E-4,1E-3,1E-4,1E-1,1E-2}, /* prn[] */
    5E-12,                      /* sclkstab */
    {3.0,0.9999,0.20},          /* thresar */
    0.0,0.0,0.05,               /* elmaskar,almaskhold,thresslip */
    30.0,30.0,30.0,             /* maxtdif,maxinno,maxgdop */
    {0},{0},{0},                /* baseline,ru,rb */
    {"",""},                    /* anttype */
    {{0}},{{0}},{0}             /* antdel,pcv,exsats */
};

    '''
    def __init__(self):
        self.mode = 2#0:single,1:dpgs,2:rtk
        self.nf = 2
        self.ionomodel = 1
        self.tropmodel = 1
        self.thresslip = 0.05
        self.maxout = 5
        self.modear = 1
        self.maxinno = 30
        self.glomodear = 1
        self.dynamics = 1
        self.prn = [ 1E-4,1E-3,1E-4,1E-1,1E-2 ,0.0]
        self.minlock = 0
        self.std = [30.0,0.03,0.3]
        self.iter = 1
        self.ena = [0,0,0,0]
        self.err = [100.0,0.003,0.003,0.0,1.0]
        self.sclkstab = 5E-12
        self.eratio = [100.0,100.0]
