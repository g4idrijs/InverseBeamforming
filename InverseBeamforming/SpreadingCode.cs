using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/// <summary>
/// Contains all the components necesary to perform simulations regaurding digital communications systems
/// and inverse beamforming
/// </summary>
namespace InverseBeamforming
{
	/// <summary>
	/// Collects many different waveform modulations into 1 class
	/// </summary>
	public partial class Modulations
	{
		/// <summary>
		/// Defines an interface to control the way spreading codes are generated and handled.
		/// </summary>
		public class SpreadingCode
		{
			/// <summary>
			/// Number of chips per symbol
			/// </summary>
			protected int _numberChips;

			/// <summary>
			/// Code matrix to use to spread the waveforms
			/// </summary>
			protected byte[,] _codeMatrix;

			/// <summary>
			/// Number of samples per symbol
			/// </summary>
			protected int _samplesPerSymbol;

			/// <summary>
			/// Number of symbols in each full length waveform
			/// </summary>
			protected int _numberSymbolsPerWaveform;

			/// <summary>
			/// Initializes the Spreading code instance. If there is no code matrix supplied, it will instantiate with gold codes of size 31
			/// </summary>
			/// <param name="numberChips">Number of chips per symbol</param>
			/// <param name="samplesPerSymbol">Number of samples per symbol</param>
			/// <param name="numberSymbolsPerWaveform">Nymber of symbols per waveform generated</param>
			/// <param name="codeMatrix">Code matrix to use when spreading (If null, uses Gold codes of size 31)</param>
			public SpreadingCode(int numberChips, int samplesPerSymbol, int numberSymbolsPerWaveform, byte[,] codeMatrix=null)
			{
				this._numberChips = numberChips;
				this._numberSymbolsPerWaveform = numberSymbolsPerWaveform;
				this._samplesPerSymbol = samplesPerSymbol;
				_codeMatrix = codeMatrix ?? getGoldCodes;
			}

			/// <summary>
			/// Gets a full waveform length spreading code, ready to multiply with the signal. 
			/// </summary>
			/// <param name="user">Index of the spreading code to use (zero based)</param>
			/// <param name="_samplesPerSymbol">Number of samples per symbol duration</param>
			/// <param name="_numChips">Number of chips per symbol duration</param>
			/// <param name="numSymbols">Number of symbols</param>
			/// <returns>Array containing the spreading code ready to be mixed with the signal</returns>
			public byte[] GetLargeSpreadingCode(int user, int numSymbols)
			{
				//Create a matrix for the spreading code
				byte[] spreadingCode = new byte[_samplesPerSymbol * numSymbols];
				int j = 0;

				//If the number of samples per chip is not an integer, throw an exception
				if (_samplesPerSymbol % _numberChips != 0)
					throw new ArgumentException("The number of samples per chip is not an integer.", "numSamples, numChips");

				//Loop through each symbol
				for (int i = 0; i < numSymbols; i++)
				{
					//Loop through each chip
					for (int k = 0; k < _numberChips; i++)
					{
						//Loop through the samples in each chip
						for (j = 0; j < _samplesPerSymbol / _numberChips; i++)
						{//              | symbol offset |       chip offset         | left in chip                Chip number
							spreadingCode[i * _samplesPerSymbol + k * _samplesPerSymbol / _numberChips + j] = _codeMatrix[user, k];
						}
					}
				}

				//Return the spreading code
				return spreadingCode;
			}

			///<summary>
			/// Applies a spreading code to a waveform
			/// </summary>
			/// <param name="waveform">Waveform to spread</param>
			/// <param name="user">Specific spreading code to use (row in the spreading code matrix)</param>
			/// <returns>Original waveform with the spreading code applied.</returns>
			public void SpreadWaveform(ref double[] waveform, int user)
			{
				//If the number of samples per chip is not an integer, throw an exception
				if (_samplesPerSymbol % _numberChips != 0)
					throw new ArgumentException("The number of samples per chip is not an integer.", "numSamples, numChips");

				//Loop through each symbol
				for (int i = 0; i < _numberSymbolsPerWaveform; i++)
				{
					//Loop through each chip
					for (int k = 0; k < _numberChips; k++)
					{
						//Loop through the samples in each chip
						for (int j = 0; j < _samplesPerSymbol / _numberChips; j++)
						{//         | symbol offset |       chip offset         | left in chip										Chip number
							waveform[i * _samplesPerSymbol + k * _samplesPerSymbol / _numberChips + j] *= ((double)_codeMatrix[user, k] * 2 - 1);
						}
					}
				}
			}

			/// <summary>
			/// Applies a spreading code to a waveform
			/// </summary>
			/// <param name="waveform">Waveform to spread</param>
			/// <param name="user">Specific spreading code to use (row in the spreading code matrix)</param>
			/// <returns>Original waveform with the spreading code applied.</returns>
			/// <remarks>In the implementation, it is functionally equivalent to SpreadWaveform. Applying the same spreading code twice gives the original signal back.</remarks>
			public void DespreadWaveform(ref double[] waveform, int user)
			{
				SpreadWaveform(ref waveform, user);
			}

			/// <summary>
			/// Generates and adds a new spread waveform with the specified power to the given waveform
			/// </summary>
			/// <param name="waveform">Waveform</param>
			/// <param name="user">User code row of the added user</param>
			/// <param name="userPower">Power of the signal of the added user</param>
			/// <param name="mod">Modulation used for the new user</param>
			public void AddSpreadWaveform(ref double[] waveform, int user, double userPower, ModulationType mod)
			{
				double[][] reference = mod.Reference;
				byte[] otherUser = mod.GenerateRandomBits();
				userPower = Math.Sqrt(userPower/mod.SignalPower);
				int index;

				//Loop through each symbol
				for (int i = 0; i < _numberSymbolsPerWaveform; i++)
				{
					//Loop through each chip
					for (int k = 0; k < _numberChips; k++)
					{
						index = i * _samplesPerSymbol + k * _samplesPerSymbol / _numberChips;
						//Loop through the samples in each chip
						for (int j = 0; j < _samplesPerSymbol / _numberChips; j++)
						{
							waveform[index + j] += userPower * reference[otherUser[i]][j]*((double)_codeMatrix[user, k] * 2 - 1);
						}
					}
				}
			}

			/// <summary>
			/// Initializes the components necesary to utilize the spreading functions
			/// </summary>
			/// <param name="codeMatrix">Code matrix to use to spread the waveforms</param>
			/// <param name="numChips">Number of chips per symbol</param>
			public void initializeSpreadingCodes(byte[,] codeMatrix, int numChips)
			{
				_numberChips = numChips;
				_codeMatrix = codeMatrix;
			}

			/// <summary>
			/// Returns the gold codes of length 31
			/// </summary>
			/// <returns>Byte array containing the gold codes of length 31</returns>
			public byte[,] getGoldCodes
			{
				get
				{
					return new byte[,] {{0,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,0,0,1},
										{1,1,0,1,1,0,1,0,0,0,0,1,1,0,0,1,0,0,1,1,1,1,1,0,1,1,1,0,0,0,1,0},
										{1,1,0,0,1,1,0,0,0,1,1,0,0,1,0,1,0,1,0,1,0,0,0,0,0,1,0,0,0,0,1,1},
										{0,0,1,0,0,0,1,0,0,1,0,0,1,1,1,0,0,0,0,1,0,0,1,1,0,1,1,0,0,1,0,0},
										{1,1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,1,0,0,1,0,1,0,1,0,0,1,0,1,0,1,1},
										{1,1,0,0,0,1,1,0,1,0,1,1,0,1,0,1,1,0,0,1,1,0,0,1,1,0,1,1,0,1,0,0},
										{0,0,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0},
										{1,1,0,1,0,1,0,1,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,0,1,1,1,1,0,1,1,1},
										{0,0,0,1,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,0,1,1,0,0,0,0,0,1,1,0,0},
										{0,0,0,1,1,0,1,0,1,1,1,0,0,0,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,0,1,1},
										{0,0,0,0,1,1,1,1,0,1,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,0,1,0,1,0,1},
										{0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,1,0,1,0,1,1,1,1,0,0,1,0,0,1},
										{1,1,1,1,0,0,1,0,1,0,0,0,0,1,1,1,1,1,1,0,0,1,0,0,0,1,1,1,0,0,0,1},
										{1,1,0,1,1,1,1,1,1,0,0,0,1,0,1,1,0,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0},
										{0,0,0,0,0,1,0,1,1,0,0,1,0,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,0,1,0},
										{0,0,1,1,0,0,0,1,1,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,0,1,1,1},
										{1,1,0,1,1,0,0,1,1,1,0,0,0,1,0,0,1,1,0,0,0,0,1,1,1,0,1,0,1,1,0,1},
										{0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,1,0,0,1,1,0,1,0,0,1,0,1,1,1,0,0,0},
										{0,0,1,0,1,0,0,0,1,0,0,1,1,1,1,0,1,1,0,1,1,0,1,0,1,0,0,1,0,0,1,1},
										{1,1,1,0,1,0,1,1,1,0,1,1,1,0,0,1,0,0,0,0,0,1,1,0,1,1,0,0,0,1,0,1},
										{1,1,1,0,1,1,0,1,1,1,1,1,0,1,1,0,1,0,1,1,1,1,1,0,0,1,1,0,1,0,0,0},
										{1,1,1,0,0,0,0,1,0,1,1,0,1,0,0,1,1,1,0,0,1,1,1,1,0,0,1,1,0,0,1,0},
										{1,1,1,1,1,0,0,0,0,1,0,1,0,1,1,1,0,0,1,0,1,1,0,1,1,0,0,0,0,1,1,0},
										{1,1,0,0,1,0,1,0,0,0,1,0,1,0,1,0,1,1,1,0,1,0,0,0,1,1,1,0,1,1,1,0},
										{0,0,1,0,1,1,1,0,1,1,0,1,0,0,0,1,0,1,1,0,0,0,1,0,0,0,1,1,1,1,1,0},
										{1,1,1,0,0,1,1,1,0,0,1,0,0,1,1,0,0,1,1,1,0,1,1,1,1,0,0,1,1,1,1,1},
										{1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,0,0,1,0,1,1,1,0,0,1,1,0,1,1,1,0,0},
										{1,1,0,1,0,0,1,1,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,1,1,0,1,0},
										{0,0,0,1,1,1,0,0,1,0,1,0,1,1,0,0,1,0,1,0,0,1,1,1,0,1,0,1,0,1,1,0},
										{0,0,0,0,0,0,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,0,1,0,1,0,0,1,1,1,1},
										{0,0,1,1,1,1,0,1,0,0,1,1,1,1,1,1,0,1,0,0,1,0,0,1,0,1,1,1,1,1,0,1},
										{1,1,0,0,0,0,0,0,1,1,1,1,1,0,1,0,0,0,1,0,0,0,0,1,0,0,0,1,1,0,0,1},
										{0,0,1,1,1,0,1,1,0,1,1,1,0,0,0,0,1,1,1,1,0,0,0,1,1,1,0,1,0,0,0,0}};
				}
			}

			/// <summary>
			/// Returns Hadamard codes of size 32
			/// </summary>
			/// <returns>Hadamard codes of length 32</returns>
			public static byte[,] getHadamardCodes
			{
				get
				{
					return new byte[,] {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
										{1,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0},
										{1,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0},
										{1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1},
										{1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0},
										{1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1},
										{1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1},
										{1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0},
										{1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0},
										{1,1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1},
										{1,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1},
										{1,1,0,0,1,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,0,1,1,0},
										{1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1},
										{1,1,0,1,0,0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,1,0,1,0},
										{1,1,1,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1,1,1,0,0},
										{1,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1},
										{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
										{1,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1},
										{1,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1},
										{1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0},
										{1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1},
										{1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1,0,1,0,1,1,0,1,0,0,1,0,1,1,0,1,0},
										{1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0},
										{1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1},
										{1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1},
										{1,1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0},
										{1,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0},
										{1,1,0,0,1,1,0,0,1,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,1,0,0,1,1,0,0,1},
										{1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0},
										{1,1,0,1,0,0,1,0,1,0,1,0,1,1,0,1,0,0,1,0,1,1,0,1,0,1,0,1,0,0,1,0,1},
										{1,1,1,0,0,0,0,1,1,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,0,0,1,1},
										{1,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0}};
				}
			}
		}
	}
}
