#include "time.h"
#include "ACMSim.h"
#define VVVF_CONTROL 1

struct InductionMachineSimulated IM;

void IM_init(){
    int i;
    for(i=0;i<5;++i){
        IM.x[i] = 0.0;
    }
    IM.rpm = 0.0;

    IM.iqs = 0.0;
    IM.ids = 0.0;

    IM.Tload = 0.0;
    IM.rpm_cmd = 0.0;
    IM.rpm_deriv_cmd = 0.0;

    IM.Lmu    = 0.4482;
    IM.Lsigma = 0.0126;

    IM.rreq   = 1.69;
    IM.rs     = 3.04;

    IM.alpha  = IM.rreq / (IM.Lmu);
    IM.Lmu_inv= 1.0/IM.Lmu;

    IM.Js = 0.0636; // Awaya92 using im.omg
    IM.npp = 2;
    IM.mu_m = IM.npp/IM.Js;

    IM.Ts  = IM_TS;

    IM.ual = 0.0;
    IM.ube = 0.0;
}

void rK5_dynamics(double t, double *x, double *fx){ //电机动态方程 [0]d电流 [1]q电流 [2]d磁链 [3]q磁链 [4]转速
    // electromagnetic model
    fx[2] = IM.rreq*x[0] - IM.alpha*x[2] - x[4]*x[3];   //psi_d'=r*id-alpha*psi_d-omega*psi_q
    fx[3] = IM.rreq*x[1] - IM.alpha*x[3] + x[4]*x[2];   //psi_q'=r*iq-alpha*psi_q-omega*psi_d
    fx[0] = (IM.ual - IM.rs*x[0] - fx[2])/IM.Lsigma;    //id'=(ual-rs*id-psi_d')/Lsigma
    fx[1] = (IM.ube - IM.rs*x[1] - fx[3])/IM.Lsigma;    //iq'=(ube-rs*iq-psi_q')/Lsigma

    // mechanical model
    IM.Tem = IM.npp*(x[1]*x[2]-x[0]*x[3]);  //Te=npp*(iq*psi_d-id*psi_q) 恒功率变换
    fx[4] = (IM.Tem - IM.Tload)*IM.mu_m;    //n'=(Te-Tl)*mu_m 电气转速
}
void rK555_Lin(double t, double *x, double hs){
    double k1[5], k2[5], k3[5], k4[5], xk[5];
    double fx[5];   //时间导数
    int i;

    rK5_dynamics(t, x, fx); // timer.t,
    for(i=0;i<5;++i){        
        k1[i] = fx[i] * hs;
        xk[i] = x[i] + k1[i]*0.5;
    }
    
    rK5_dynamics(t, xk, fx); // timer.t+hs/2., 
    for(i=0;i<5;++i){        
        k2[i] = fx[i] * hs;
        xk[i] = x[i] + k2[i]*0.5;
    }
    
    rK5_dynamics(t, xk, fx); // timer.t+hs/2., 
    for(i=0;i<5;++i){        
        k3[i] = fx[i] * hs;
        xk[i] = x[i] + k3[i];
    }
    
    rK5_dynamics(t, xk, fx); // timer.t+hs, 
    for(i=0;i<5;++i){        
        k4[i] = fx[i] * hs;
        x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6.0;
    }
}
int machine_simulation(){
    rK555_Lin(CTRL.timebase, IM.x, IM.Ts);  //龙科库塔法 5个变量 线性  时基、电机状态、仿真步长

    IM.ids = IM.x[0];   //id电流
    IM.iqs = IM.x[1];   //iq电流
    IM.rpm = IM.x[4] * 60 / (2 * M_PI * IM.npp);    //x[4]we(rad/s)->(rpm)   we=npp*wm=2*pi*n/60

    if(isNumber(IM.rpm))    //转速不跑飞
        return false;
    else
        return true;
}
void measurement(){
    US_C(0) = CTRL.ual; //观测器电压就是控制器输出电压 电压给定就是测量电压 当前电压就是上一步电压
    US_C(1) = CTRL.ube;
    US_P(0) = US_C(0);
    US_P(1) = US_C(1);

    IS_C(0) = IM.ids;
    IS_C(1) = IM.iqs;

    im.omg = IM.x[4];
}
void inverter_model(){
    IM.ual = CTRL.ual;  //电机电压就是控制器电压
    IM.ube = CTRL.ube;
}
int main(){
    printf("NUMBER_OF_LINES: %d\n\n", NUMBER_OF_LINES); //仿真步长

    /* Initialization */
    IM_init();  //电机初始化
    CTRL_init();    //控制器初始化

    FILE *fw;   //写入文件，指针
    fw = fopen("algorithm.dat", "w");   //写入文件

    /* MAIN LOOP */
    clock_t  begin, end;
    begin = clock();    //运行时间
    int _; // _ for the outer iteration
    int dfe=0; // dfe for down frequency execution
    for(_=0;_<NUMBER_OF_LINES;++_){

        /* Command and Load Torque */
        IM.rpm_cmd = 50; // rpm 转速指令
        IM.Tload = 10; // Nm 负载指令

        /* Simulated IM */
        if(machine_simulation()){ 
            printf("Break the loop.\n");
            break;
        } 

        if(++dfe==DOWN_FREQ_EXE){   //dsp里面 采样 控制器执行频率低于仿真频率，(e.g.仿真8kHz，控制器4kHz)
            dfe = 0;

            /* Time */
            CTRL.timebase += TS;

            measurement();  //电流采样

            // observation();

            write_data_to_file(fw);

            #if VVVF_CONTROL == TRUE
                #define VF_RATIO 18 //18.0 // 8 ~ 18 shows saturated phenomenon
                double freq = 2; // 0.15 ~ 0.5 ~ 2 （0.1时电压李萨茹就变成一个圆了）
                double volt = VF_RATIO*freq;
                CTRL.ual = volt*cos(2*M_PI*freq*CTRL.timebase);
                CTRL.ube = volt*sin(2*M_PI*freq*CTRL.timebase);
            #else
                control();
            #endif
        }

        inverter_model();   //控制器给定就是逆变器输出
    }
    end = clock(); printf("The simulation in C costs %g sec.\n", (double)(end - begin)/CLOCKS_PER_SEC);
    fclose(fw);

    /* Fade out */
    system("python ./ACMPlot.py"); 
    // getch();
    // system("pause");
    // system("exit");
    return 0; 
}


/* Utility */
void write_data_to_file(FILE *fw){
    static int j=0,jj=0; // j,jj for down sampling

    // if(CTRL.timebase>20)
    {
    if(++j == 10)   //控制器每10步输出1步数据，写入
    {
        j=0;
        fprintf(fw, "%g,%g,%g,%g,%g\n",
                IM.x[0], IM.x[1], IM.x[2], IM.x[3], IM.x[4]/IM.npp/2/M_PI*60 // from elec.rad/s to rpm
                );
    }
    }
}

bool isNumber(double x){
    // This looks like it should always be true, 
    // but it's false if x is a NaN (1.#QNAN0).
    return (x == x); 
    // see https://www.johndcook.com/blog/IEEE_exceptions_in_cpp/ cb: https://stackoverflow.com/questions/347920/what-do-1-inf00-1-ind00-and-1-ind-mean
}



