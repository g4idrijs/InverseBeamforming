using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming
{
	public partial class Filter
	{
		/// <summary>
		/// Gets the filter coefficients for the DS filter
		/// </summary>
		/// <param name="IIRCoeff">Structure to fill wil IIR coeficients for the filter</param>
		public static void getIIRCoefficients_Gold_DS(out IIR_Filter.TIIRCoeff IIRCoeff)
		{
			IIRCoeff = new IIR_Filter.TIIRCoeff(true);
			IIRCoeff.NumSections = 8;
			IIRCoeff.b0[0] = 1.162221116541298e-12;
			IIRCoeff.b1[0] = 2.3473882338099569e-12;
			IIRCoeff.b2[0] = 1.1853000077051452e-12;
			IIRCoeff.a0[0] = 1;
			IIRCoeff.a1[0] = -1.6698316775306086;
			IIRCoeff.a2[0] = 0.9246205671829052;
			IIRCoeff.b0[1] = 1;
			IIRCoeff.b1[1] = 2.0080755879068581;
			IIRCoeff.b2[1] = 1.0081893423258401;
			IIRCoeff.a0[1] = 1;
			IIRCoeff.a1[1] = -1.7908983344318448;
			IIRCoeff.a2[1] = 1.0772534421510178;
			IIRCoeff.b0[2] = 1;
			IIRCoeff.b1[2] = 1.9917793501816372;
			IIRCoeff.b2[2] = 0.9918922686263234;
			IIRCoeff.a0[2] = 1;
			IIRCoeff.a1[2] = -1.8102207227019464;
			IIRCoeff.a2[2] = 1.0649079622780933;
			IIRCoeff.b0[3] = 1;
			IIRCoeff.b1[3] = 1.9804018302658404;
			IIRCoeff.b2[3] = 0.98051415444574963;
			IIRCoeff.a0[3] = 1;
			IIRCoeff.a1[3] = -1.6567402382608996;
			IIRCoeff.a2[3] = 0.94080669142659423;
			IIRCoeff.b0[4] = 1;
			IIRCoeff.b1[4] = -2.0200906079675218;
			IIRCoeff.b2[4] = 1.0202089824983394;
			IIRCoeff.a0[4] = 1;
			IIRCoeff.a1[4] = -1.7196489021938883;
			IIRCoeff.a2[4] = 0.9518172540106159;
			IIRCoeff.b0[5] = 1;
			IIRCoeff.b1[5] = -1.9800339583548481;
			IIRCoeff.b2[5] = 0.98015056882639628;
			IIRCoeff.a0[5] = 1;
			IIRCoeff.a1[5] = -1.7379145908672173;
			IIRCoeff.a2[5] = 1.0442075461657843;
			IIRCoeff.b0[6] = 1;
			IIRCoeff.b1[6] = -2.0082333132048422;
			IIRCoeff.b2[6] = 1.0083511654118005;
			IIRCoeff.a0[6] = 1;
			IIRCoeff.a1[6] = -1.6839887555861679;
			IIRCoeff.a2[6] = 0.98871758363282636;
			IIRCoeff.b0[7] = 1;
			IIRCoeff.b1[7] = -1.991642120472791;
			IIRCoeff.b2[7] = 0.9917592417354667;
			IIRCoeff.a0[7] = 1;
			IIRCoeff.a1[7] = -1.7798530468574201;
			IIRCoeff.a2[7] = 1.0111722814023301;
		}

		/// <summary>
		///Gets the numerator coefficients for the DS filter.
		/// </summary>
		/// <param name="numerator">Array to hold the coefficients of the numerator</param>
		public static void getNumeratorFilterCoefficients_Gold_DS(out double[] numerator)
		{
			numerator = new double[]{1.162221116541298e-12,0,-9.2977689323303842e-12,0,3.2542191263156341e-11,0,-6.5084382526312683e-11,0,8.1355478157890857e-11,0,-6.5084382526312683e-11,0,3.2542191263156341e-11,0,-9.2977689323303842e-12,0,1.162221116541298e-12};
		}

		/// <summary>
		///Gets the denominator coefficients for the DS filter.
		/// </summary>
		/// <param name="denominator">Array to hold the coefficients of the numerator</param>
		public static void getDenominatorFilterCoefficients_Gold_DS(out double[] denominator)
		{
			denominator = new double[]{1,-13.849096268429992,91.902948840949989,-387.3655983184558,1159.5495248195007,-2611.9589880711501,4577.4811682419004,-6364.2461801519275,7092.9861476309752,-6357.5302657991833,4567.8254288960943,-2603.6988479751822,1154.6627765037933,-385.32605635797233,91.322594494097814,-13.747119085342884,0.99158905960120747};
		}

		/// <summary>
		/// Gets the filter coefficients for the RF filter
		/// </summary>
		/// <param name="IIRCoeff">Structure to fill wil IIR coeficients for the filter</param>
		public static void getIIRCoefficients_Gold_RF(out IIR_Filter.TIIRCoeff IIRCoeff)
		{
			IIRCoeff = new IIR_Filter.TIIRCoeff(true);
			IIRCoeff.NumSections = 8;
			IIRCoeff.b0[0] = -4.9370624554234931e-12;
			IIRCoeff.b1[0] = -9.8741800442909721e-12;
			IIRCoeff.b2[0] = -4.9367935558984279e-12;
			IIRCoeff.a0[0] = 1;
			IIRCoeff.a1[0] = -1.6402378402218392;
			IIRCoeff.a2[0] = 0.89726802178109744;
			IIRCoeff.b0[1] = 1;
			IIRCoeff.b1[1] = -2.0001216849367198;
			IIRCoeff.b2[1] = 1.0000303639368076;
			IIRCoeff.a0[1] = 1;
			IIRCoeff.a1[1] = -1.6859180998887944;
			IIRCoeff.a2[1] = 0.91538781709900319;
			IIRCoeff.b0[2] = 1;
			IIRCoeff.b1[2] = 2.0114606189783757;
			IIRCoeff.b2[2] = 1.0115263358648694;
			IIRCoeff.a0[2] = 1;
			IIRCoeff.a1[2] = -1.6325059818643695;
			IIRCoeff.a2[2] = 0.91977534811672035;
			IIRCoeff.b0[3] = 1;
			IIRCoeff.b1[3] = 1.9999888169901794;
			IIRCoeff.b2[3] = 1.0000544898205288;
			IIRCoeff.a0[3] = 1;
			IIRCoeff.a1[3] = -1.7650254610509331;
			IIRCoeff.a2[3] = 1.0400584514823925;
			IIRCoeff.b0[4] = 1;
			IIRCoeff.b1[4] = 1.9885393967747151;
			IIRCoeff.b2[4] = 0.98860498553617171;
			IIRCoeff.a0[4] = 1;
			IIRCoeff.a1[4] = -1.6620429876074099;
			IIRCoeff.a2[4] = 0.96846319319726093;
			IIRCoeff.b0[5] = 1;
			IIRCoeff.b1[5] = -2.0135138317309247;
			IIRCoeff.b2[5] = 1.0136059759089591;
			IIRCoeff.a0[5] = 1;
			IIRCoeff.a1[5] = -1.7494355296848343;
			IIRCoeff.a2[5] = 0.97022394029812009;
			IIRCoeff.b0[6] = 1;
			IIRCoeff.b1[6] = -1.9864861504219293;
			IIRCoeff.b2[6] = 0.98657665029542208;
			IIRCoeff.a0[6] = 1;
			IIRCoeff.a1[6] = -1.785487411362874;
			IIRCoeff.a2[6] = 1.0249157178723367;
			IIRCoeff.b0[7] = 1;
			IIRCoeff.b1[7] = -1.9998783329104248;
			IIRCoeff.b2[7] = 0.99996964855309656;
			IIRCoeff.a0[7] = 1;
			IIRCoeff.a1[7] = -1.713441054504441;
			IIRCoeff.a2[7] = 1.0170853491624059;
		}

		/// <summary>
		///Gets the numerator coefficients for the RF filter.
		/// </summary>
		/// <param name="numerator">Array to hold the coefficients of the numerator</param>
		public static void getNumeratorFilterCoefficients_Gold_RF(out double[] numerator)
		{
			numerator = new double[]{-4.9370624554234931e-12,-0,3.9496499643387944e-11,-0,-1.3823774875185783e-10,-0,2.7647549750371566e-10,-0,-3.4559437187964451e-10,-0,2.7647549750371566e-10,-0,-1.3823774875185783e-10,-0,3.9496499643387944e-11,-0,-4.9370624554234931e-12};
		}

		/// <summary>
		///Gets the denominator coefficients for the RF filter.
		/// </summary>
		/// <param name="denominator">Array to hold the coefficients of the numerator</param>
		public static void getDenominatorFilterCoefficients_Gold_RF(out double[] denominator)
		{
			denominator = new double[]{1,-13.634094366185479,89.067518542242482,-369.55646143364493,1088.9470330959643,-2414.5236862988809,4165.150964718775,-5700.1010880057238,6253.0079749799352,-5516.525694846131,3901.1885043240381,-2188.6725587436249,955.29889531493677,-313.75931580824039,73.184428232242283,-10.842001092097329,0.76960430657535384};
		}

		/// <summary>
		/// Gets the filter coefficients for the DS filter
		/// </summary>
		/// <param name="IIRCoeff">Structure to fill wil IIR coeficients for the filter</param>
		public static void getIIRCoefficients_Hadamard_DS(out IIR_Filter.TIIRCoeff IIRCoeff)
		{
			IIRCoeff = new IIR_Filter.TIIRCoeff(true);
			IIRCoeff.NumSections = 8;
			IIRCoeff.b0[0] = -4.0651683029553006e-13;
			IIRCoeff.b1[0] = -8.1310255659040189e-13;
			IIRCoeff.b2[0] = -4.0655387934558122e-13;
			IIRCoeff.a0[0] = 1;
			IIRCoeff.a1[0] = -1.7832183120297342;
			IIRCoeff.a2[0] = 1.0631487775929527;
			IIRCoeff.b0[1] = 1;
			IIRCoeff.b1[1] = 2.0125126630222994;
			IIRCoeff.b2[1] = 1.0125920468714575;
			IIRCoeff.a0[1] = 1;
			IIRCoeff.a1[1] = -1.6777119102933238;
			IIRCoeff.a2[1] = 0.93911188739822715;
			IIRCoeff.b0[2] = 1;
			IIRCoeff.b1[2] = 1.9998304827674178;
			IIRCoeff.b2[2] = 0.99990877307922621;
			IIRCoeff.a0[2] = 1;
			IIRCoeff.a1[2] = -1.7953792018617887;
			IIRCoeff.a2[2] = 1.048866421777481;
			IIRCoeff.b0[3] = 1;
			IIRCoeff.b1[3] = 1.9874873753739775;
			IIRCoeff.b2[3] = 0.9875646372332092;
			IIRCoeff.a0[3] = 1;
			IIRCoeff.a1[3] = -1.7130928237405842;
			IIRCoeff.a2[3] = 0.95383218729524777;
			IIRCoeff.b0[4] = 1;
			IIRCoeff.b1[4] = -2.0159765017114117;
			IIRCoeff.b2[4] = 1.0160511284795004;
			IIRCoeff.a0[4] = 1;
			IIRCoeff.a1[4] = -1.6741512232361586;
			IIRCoeff.a2[4] = 0.95794894095238581;
			IIRCoeff.b0[5] = 1;
			IIRCoeff.b1[5] = -1.9839024496907354;
			IIRCoeff.b2[5] = 0.98397845027984865;
			IIRCoeff.a0[5] = 1;
			IIRCoeff.a1[5] = -1.7419527894856679;
			IIRCoeff.a2[5] = 1.0390803871872047;
			IIRCoeff.b0[6] = 1;
			IIRCoeff.b1[6] = -2.0067101821484483;
			IIRCoeff.b2[6] = 1.0067852917139737;
			IIRCoeff.a0[6] = 1;
			IIRCoeff.a1[6] = -1.6991717298485969;
			IIRCoeff.a2[6] = 0.9968524223369527;
			IIRCoeff.b0[7] = 1;
			IIRCoeff.b1[7] = -1.9934108664493957;
			IIRCoeff.b2[7] = 0.99348654100582034;
			IIRCoeff.a0[7] = 1;
			IIRCoeff.a1[7] = -1.7644696449194333;
			IIRCoeff.a2[7] = 1.0005334285331222;
		}

		/// <summary>
		///Gets the numerator coefficients for the DS filter.
		/// </summary>
		/// <param name="numerator">Array to hold the coefficients of the numerator</param>
		public static void getNumeratorFilterCoefficients_Hadamard_DS(out double[] numerator)
		{
			numerator = new double[]{-4.0651683029553006e-13,0,3.2521346423642405e-12,0,-1.1382471248274839e-11,0,2.2764942496549679e-11,0,-2.8456178120687101e-11,0,2.2764942496549679e-11,0,-1.1382471248274839e-11,0,3.2521346423642405e-12,0,-4.0651683029553006e-13};
		}

		/// <summary>
		///Gets the denominator coefficients for the DS filter.
		/// </summary>
		/// <param name="denominator">Array to hold the coefficients of the numerator</param>
		public static void getDenominatorFilterCoefficients_Hadamard_DS(out double[] denominator)
		{
			denominator = new double[]{1,-13.849147635415274,91.903630647864645,-387.3699093146588,1159.5667321529022,-2612.0074416858051,4577.5830720829672,-6364.4114823540403,7093.19670579993,-6357.742591330888,4567.9949400788037,-2603.8051372771574,1154.7141997935737,-385.34464779581492,91.327339822607115,-13.747884476189382,0.99164795138241879};
		}

		/// <summary>
		/// Gets the filter coefficients for the RF filter
		/// </summary>
		/// <param name="IIRCoeff">Structure to fill wil IIR coeficients for the filter</param>
		public static void getIIRCoefficients_Hadamard_RF(out IIR_Filter.TIIRCoeff IIRCoeff)
		{
			IIRCoeff = new IIR_Filter.TIIRCoeff(true);
			IIRCoeff.NumSections = 8;
			IIRCoeff.b0[0] = -6.8352739083232691e-13;
			IIRCoeff.b1[0] = -1.3795729091285681e-12;
			IIRCoeff.b2[0] = -6.9611284206494148e-13;
			IIRCoeff.a0[0] = 1;
			IIRCoeff.a1[0] = -1.6459050566742155;
			IIRCoeff.a2[0] = 0.90241581968705176;
			IIRCoeff.b0[1] = 1;
			IIRCoeff.b1[1] = 2.0074228424265268;
			IIRCoeff.b2[1] = 1.0075204676349672;
			IIRCoeff.a0[1] = 1;
			IIRCoeff.a1[1] = -1.6890686536519219;
			IIRCoeff.a2[1] = 0.92010046662978784;
			IIRCoeff.b0[2] = 1;
			IIRCoeff.b1[2] = 1.9923465855354747;
			IIRCoeff.b2[2] = 0.99244298170369083;
			IIRCoeff.a0[2] = 1;
			IIRCoeff.a1[2] = -1.637328923940363;
			IIRCoeff.a2[2] = 0.92211733322310696;
			IIRCoeff.b0[3] = 1;
			IIRCoeff.b1[3] = 1.9819165620483412;
			IIRCoeff.b2[3] = 0.98201208957084529;
			IIRCoeff.a0[3] = 1;
			IIRCoeff.a1[3] = -1.6630708997040147;
			IIRCoeff.a2[3] = 0.96640529559442589;
			IIRCoeff.b0[4] = 1;
			IIRCoeff.b1[4] = -2.0187175863977771;
			IIRCoeff.b2[4] = 1.0188204206380649;
			IIRCoeff.a0[4] = 1;
			IIRCoeff.a1[4] = -1.7567422542244244;
			IIRCoeff.a2[4] = 1.0314950526125641;
			IIRCoeff.b0[5] = 1;
			IIRCoeff.b1[5] = -1.9814750855031913;
			IIRCoeff.b2[5] = 0.98157538282580403;
			IIRCoeff.a0[5] = 1;
			IIRCoeff.a1[5] = -1.7480922187809054;
			IIRCoeff.a2[5] = 0.97132527374718669;
			IIRCoeff.b0[6] = 1;
			IIRCoeff.b1[6] = -2.0076167890276588;
			IIRCoeff.b2[6] = 1.0077188796252752;
			IIRCoeff.a0[6] = 1;
			IIRCoeff.a1[6] = -1.7795645299715708;
			IIRCoeff.a2[6] = 1.0203738845556298;
			IIRCoeff.b0[7] = 1;
			IIRCoeff.b1[7] = -1.9921905390713719;
			IIRCoeff.b2[7] = 0.9922915787372335;
			IIRCoeff.a0[7] = 1;
			IIRCoeff.a1[7] = -1.7088762376548443;
			IIRCoeff.a2[7] = 1.0107522870505476;
		}

		/// <summary>
		///Gets the numerator coefficients for the RF filter.
		/// </summary>
		/// <param name="numerator">Array to hold the coefficients of the numerator</param>
		public static void getNumeratorFilterCoefficients_Hadamard_RF(out double[] numerator)
		{
			numerator = new double[]{-6.8352739083232691e-13,0,5.4682191266586153e-12,0,-1.9138766943305151e-11,0,3.8277533886610302e-11,0,-4.7846917358262876e-11,0,3.8277533886610302e-11,0,-1.9138766943305151e-11,0,5.4682191266586153e-12,0,-6.8352739083232691e-13};
		}

		/// <summary>
		///Gets the denominator coefficients for the RF filter.
		/// </summary>
		/// <param name="denominator">Array to hold the coefficients of the numerator</param>
		public static void getDenominatorFilterCoefficients_Hadamard_RF(out double[] denominator)
		{
			denominator = new double[]{1,-13.628648774602238,88.996176693812231,-369.11144803288579,1087.1953128776229,-2409.6605939391575,4155.0696257077561,-5683.9844493445471,6232.7787328803224,-5496.4264864242232,3885.3791516741344,-2178.9063661702485,950.64402717538985,-312.10135156166245,72.767503055352847,-10.77574371337945,0.76458072201601168};
		}
	}
}
