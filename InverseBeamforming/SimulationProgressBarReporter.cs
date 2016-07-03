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
		public class SimulationProgressBarReporter : IObserver<SingleSimulation>
		{
			public event SimulationProgressUpdated SimulationProgressBarUpdatedEvent;
			public delegate void SimulationProgressUpdated(SimulationProgressBarReporter s, EventArgs e);

			public event SimulationProgressUpdated SimulationErrorOccuredEvent;
			public delegate void SimulationErrorOccured(SimulationProgressBarReporter s, EventArgs e);

			/// <summary>
			/// Percentage of the simulation that is completed
			/// </summary>
			public double Progress
			{
				get { return this._progress; }
			}
			private double _progress;

			/// <summary>
			/// Length of time the simulation has been running, or if it is done, the length of time it was running.
			/// </summary>
			public TimeSpan LengthOfSimulation { get { return this._lengthOfSimulation; } }
			private TimeSpan _lengthOfSimulation;

			/// <summary>
			/// Information to unsubscribe from the simulation
			/// </summary>
			private IDisposable unsubscriber;

			/// <summary>
			/// Filename of the log file to write the intermediate results to
			/// </summary>
			private string _instName;

			/// <summary>
			/// Start and end times of the simulation
			/// </summary>
			private DateTime startTime, endTime;
			
			/// <summary>
			/// Construct a new instances of the SimulationReporter class
			/// </summary>
			/// <param name="simName">Name of the logfile to use for this reporter</param>
			public SimulationProgressBarReporter(string simName)
			{
				this._instName = simName;
				this._lengthOfSimulation = new TimeSpan(0, 0, 0);
			}

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
			/// Sets up the file to log the simulation
			/// </summary>
			public virtual void OnStart()
			{
				//Get the starting time
				startTime = DateTime.Now;

				//Reset the progress bar to 0% complete
				this.setProgressBar(0);

				//Tell subscribers that the progress bar was updated
				SimulationProgressBarUpdatedEvent(this, new EventArgs());
			}

			/// <summary>
			/// Report the status of the simulation to the logfile
			/// </summary>
			/// <param name="sim">Simulation to report the status of</param>
			public virtual void OnNext(SingleSimulation sim)
			{
				//Get the intermediate results from the simulation
				IntermediateSimResults isr = sim.isr;

				//Update the progress bar
				this.setProgressBar(isr.PercentErrorHad);
				this._lengthOfSimulation = DateTime.Now - startTime;

				//Tell subscribers that the progress bar was updated
				SimulationProgressBarUpdatedEvent(this, new EventArgs());
			}

			/// <summary>
			/// Report an error in the simulation to the logfile
			/// </summary>
			/// <param name="error">Exception describing the error</param>
			public virtual void OnError(Exception error)
			{
				//Tell subscribers that an error occured
				SimulationErrorOccuredEvent(this, new ErrorOccuredEventArgs(error));
			}

			/// <summary>
			/// Write completion information to the logfile and unsuscribe from the simulation
			/// </summary>
			public virtual void OnCompleted()
			{
				//Get the end time of the simulation
				endTime = DateTime.Now;

				this.setProgressBar(100);
				this._lengthOfSimulation = endTime - startTime;

				//Tell subscribers that the progress bar was updated
				SimulationProgressBarUpdatedEvent(this, new EventArgs());

				//Unsubscribe from the simulation
				this.Unsubscribe();
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
			/// Private setter for the progress bar
			/// </summary>
			/// <param name="value">New value for the bar</param>
			private void setProgressBar(double value)
			{
				if (value >= 0 && value <= 100)
					this._progress = value;
				else
					throw new ArgumentOutOfRangeException("Progress", "Progress must be a percentage between 0 and 100 inclusive.");
			}
		}
	}

	/// <summary>
	/// Event args class that can pass an exception
	/// </summary>
	public class ErrorOccuredEventArgs : EventArgs
	{
		public Exception e;

		public ErrorOccuredEventArgs(Exception e) : base()
		{
			this.e = e;
		}
	}
}
