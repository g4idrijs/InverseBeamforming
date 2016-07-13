using System;
using System.Collections.Generic;


namespace InverseBeamforming
{
	public partial class Filter
	{
		public static void getIIR_Coefficients_Gold_DS(out IIR_Filter.TIIRCoeff coefs)
		{
			coefs = new IIR_Filter.TIIRCoeff(true);
			coefs.a0[1] = 0.00070413321580150461;
			coefs.a0[2] = 0.00070413321580150461;
			coefs.a0[3] = 0.00070393452932227517;
			coefs.a0[4] = 0.00070393452932227517;
			coefs.a0[5] = 0.00070379713244000128;
			coefs.a0[6] = 0.00070379713244000128;
			coefs.a0[7] = 0.00070374811553745875;
			coefs.a1[1] = -1.7310923922594577;
			coefs.a1[2] = -1.7324659955166632;
			coefs.a1[3] = -1.7307395977902;
			coefs.a1[4] = -1.7318417632443899;
			coefs.a1[5] = -1.7306469944097309;
			coefs.a1[6] = -1.7312588884159359;
			coefs.a1[7] = -1.7308324675544136;
			coefs.a2[1] = 0.9996862587110722;
			coefs.a2[2] = 0.99968700391354992;
			coefs.a2[3] = 0.99912137137990831;
			coefs.a2[4] = 0.99912304488606463;
			coefs.a2[5] = 0.9987311306864024;
			coefs.a2[6] = 0.99873247220986539;
			coefs.a2[7] = 0.99859250376892517;
			coefs.b0[1] = 1;
			coefs.b0[2] = 1;
			coefs.b0[3] = 1;
			coefs.b0[4] = 1;
			coefs.b0[5] = 1;
			coefs.b0[6] = 1;
			coefs.b0[7] = 1;
			coefs.b1[1] = 0;
			coefs.b1[2] = 0;
			coefs.b1[3] = 0;
			coefs.b1[4] = 0;
			coefs.b1[5] = 0;
			coefs.b1[6] = 0;
			coefs.b1[7] = 0;
			coefs.b2[1] = -1;
			coefs.b2[2] = -1;
			coefs.b2[3] = -1;
			coefs.b2[4] = -1;
			coefs.b2[5] = -1;
			coefs.b2[6] = -1;
			coefs.b2[7] = -1;
		}
		public static void getIIR_Coefficients_Gold_RF(out IIR_Filter.TIIRCoeff coefs)
		{
			coefs = new IIR_Filter.TIIRCoeff(true);
			coefs.a0[1] = 0.021725833739074118;
			coefs.a0[2] = 0.021725833739074118;
			coefs.a0[3] = 0.021539942426096515;
			coefs.a0[4] = 0.021539942426096515;
			coefs.a0[5] = 0.021414060129567619;
			coefs.a0[6] = 0.021414060129567619;
			coefs.a0[7] = 0.021369664771290865;
			coefs.a1[1] = -1.7019811484262308;
			coefs.a1[2] = -1.7449473174581205;
			coefs.a1[3] = -1.6913419189584451;
			coefs.a1[4] = -1.7263759324434287;
			coefs.a1[5] = -1.6892102394364403;
			coefs.a1[6] = -1.7088692439197495;
			coefs.a1[7] = -1.6955905028879716;
			coefs.a2[1] = 0.98997846909011589;
			coefs.a2[2] = 0.99068822751544738;
			coefs.a2[3] = 0.97236003479660005;
			coefs.a2[4] = 0.97392868802171106;
			coefs.a2[5] = 0.96079297678572728;
			coefs.a2[6] = 0.96203728461035953;
			coefs.a2[7] = 0.95726067045741825;
			coefs.b0[1] = 1;
			coefs.b0[2] = 1;
			coefs.b0[3] = 1;
			coefs.b0[4] = 1;
			coefs.b0[5] = 1;
			coefs.b0[6] = 1;
			coefs.b0[7] = 1;
			coefs.b1[1] = 0;
			coefs.b1[2] = 0;
			coefs.b1[3] = 0;
			coefs.b1[4] = 0;
			coefs.b1[5] = 0;
			coefs.b1[6] = 0;
			coefs.b1[7] = 0;
			coefs.b2[1] = -1;
			coefs.b2[2] = -1;
			coefs.b2[3] = -1;
			coefs.b2[4] = -1;
			coefs.b2[5] = -1;
			coefs.b2[6] = -1;
			coefs.b2[7] = -1;
		}
		public static void getIIR_Coefficients_Hadamard_DS(out IIR_Filter.TIIRCoeff coefs)
		{
			coefs = new IIR_Filter.TIIRCoeff(true);
			coefs.a0[1] = 0.00069918304950063398;
			coefs.a0[2] = 0.00069918304950063398;
			coefs.a0[3] = 0.00069898714599631474;
			coefs.a0[4] = 0.00069898714599631474;
			coefs.a0[5] = 0.00069885167293218199;
			coefs.a0[6] = 0.00069885167293218199;
			coefs.a0[7] = 0.00069880334222482043;
			coefs.a1[1] = -1.7310991339141522;
			coefs.a1[2] = -1.7324630760212818;
			coefs.a1[3] = -1.7307488178286681;
			coefs.a1[4] = -1.731843226949493;
			coefs.a1[5] = -1.7306568598311589;
			coefs.a1[6] = -1.7312644460704316;
			coefs.a1[7] = -1.7308410239235303;
			coefs.a2[1] = 0.99968846696256808;
			coefs.a2[2] = 0.99968920172147402;
			coefs.a2[3] = 0.99912755238269946;
			coefs.a2[4] = 0.99912920244988612;
			coefs.a2[5] = 0.99874005148952527;
			coefs.a2[6] = 0.99874137422446563;
			coefs.a2[7] = 0.99860239331555034;
			coefs.b0[1] = 1;
			coefs.b0[2] = 1;
			coefs.b0[3] = 1;
			coefs.b0[4] = 1;
			coefs.b0[5] = 1;
			coefs.b0[6] = 1;
			coefs.b0[7] = 1;
			coefs.b1[1] = 0;
			coefs.b1[2] = 0;
			coefs.b1[3] = 0;
			coefs.b1[4] = 0;
			coefs.b1[5] = 0;
			coefs.b1[6] = 0;
			coefs.b1[7] = 0;
			coefs.b2[1] = -1;
			coefs.b2[2] = -1;
			coefs.b2[3] = -1;
			coefs.b2[4] = -1;
			coefs.b2[5] = -1;
			coefs.b2[6] = -1;
			coefs.b2[7] = -1;
		}
		public static void getIIR_Coefficients_Hadamard_RF(out IIR_Filter.TIIRCoeff coefs)
		{
			coefs = new IIR_Filter.TIIRCoeff(true);
			coefs.a0[1] = 0.022266280124460171;
			coefs.a0[2] = 0.022266280124460171;
			coefs.a0[3] = 0.022071113020399803;
			coefs.a0[4] = 0.022071113020399803;
			coefs.a0[5] = 0.021939019125336238;
			coefs.a0[6] = 0.021939019125336238;
			coefs.a0[7] = 0.021892446462511397;
			coefs.a1[1] = -1.7012202640216549;
			coefs.a1[2] = -1.7452703745995208;
			coefs.a1[3] = -1.6903230217699421;
			coefs.a1[4] = -1.7262550793803035;
			coefs.a1[5] = -1.6881569639409681;
			coefs.a1[6] = -1.7083250317154333;
			coefs.a1[7] = -1.6947127137144216;
			coefs.a2[1] = 0.98972023876610726;
			coefs.a2[2] = 0.99046576476670678;
			coefs.a2[3] = 0.97165880295864016;
			coefs.a2[4] = 0.97330586662077301;
			coefs.a2[5] = 0.95981624579774649;
			coefs.a2[6] = 0.96112243501656613;
			coefs.a2[7] = 0.95621510707497714;
			coefs.b0[1] = 1;
			coefs.b0[2] = 1;
			coefs.b0[3] = 1;
			coefs.b0[4] = 1;
			coefs.b0[5] = 1;
			coefs.b0[6] = 1;
			coefs.b0[7] = 1;
			coefs.b1[1] = 0;
			coefs.b1[2] = 0;
			coefs.b1[3] = 0;
			coefs.b1[4] = 0;
			coefs.b1[5] = 0;
			coefs.b1[6] = 0;
			coefs.b1[7] = 0;
			coefs.b2[1] = -1;
			coefs.b2[2] = -1;
			coefs.b2[3] = -1;
			coefs.b2[4] = -1;
			coefs.b2[5] = -1;
			coefs.b2[6] = -1;
			coefs.b2[7] = -1;
		}
	}
}
