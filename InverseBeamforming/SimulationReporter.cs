using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming
{
	/// <summary>
	/// Holds methods that create and report on simulations of Bit error rate
	/// </summary>
	public partial class Simulations
	{
		/// <summary>
		/// Class that describes the process of being subscribed to a simulation
		/// </summary>
		public abstract class SimulationReporter : IObserver<SingleSimulation>
		{
			/// <summary>
			/// Information to unsubscribe from the simulation
			/// </summary>
			protected IDisposable unsubscriber;

			/// <summary>
			/// Filename of the log file to write the intermediate results to
			/// </summary>
			protected string _logFilename;

			/// <summary>
			/// Start and end times of the simulation
			/// </summary>
			protected DateTime startTime, endTime;
			
			/// <summary>
			/// Subscribe to the simulation
			/// </summary>
			/// <param name="provider">Requesting observer of the simulation</param>
			public virtual void Subscribe(IObservable<SingleSimulation> provider)
			{
				if (provider != null)
					unsubscriber = provider.Subscribe(this);
			}

			/// <summary>
			/// Unsubscribe from the simulation
			/// </summary>
			public virtual void Unsubscribe()
			{
				if (unsubscriber != null)
					unsubscriber.Dispose();
			}

			/// <summary>
			/// Sets up the file to log the simulation
			/// </summary>
			public abstract void OnStart();

			/// <summary>
			/// Report the status of the simulation to the logfile
			/// </summary>
			/// <param name="sim">Simulation to report the status of</param>
			public abstract void OnNext(SingleSimulation sim);

			/// <summary>
			/// Report an error in the simulation to the logfile
			/// </summary>
			/// <param name="error">Exception describing the error</param>
			public abstract void OnError(Exception error);

			/// <summary>
			/// Write completion information to the logfile and unsuscribe from the simulation
			/// </summary>
			public abstract void OnCompleted();
		}
	}
}
