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
	/// Holds methods that create and report on simulations of Bit error rate
	/// </summary>
	public partial class Simulations
	{
		/// <summary>
		/// Class that provides methods to track a simulation
		/// </summary>
		public class SimulationTracker : IObservable<SingleSimulation>
		{
			/// <summary>
			/// List of observers of the simulation
			/// </summary>
			protected List<IObserver<SingleSimulation>> observers;

			/// <summary>
			/// Construct a new tracker
			/// </summary>
			public SimulationTracker()
			{
				observers = new List<IObserver<SingleSimulation>>();
			}

			/// <summary>
			/// Provides a method to subscribe to the simulation
			/// </summary>
			/// <param name="observer">Requesting observer of the simulation</param>
			/// <returns>IDisposable object that can unsubscribe from the simulation</returns>
			public IDisposable Subscribe(IObserver<SingleSimulation> observer)
			{
				if (!observers.Contains(observer))
					observers.Add(observer);

				return new Unsubsriber(observers, observer);
			}

			/// <summary>
			/// Class details how to unsubscribe from the simulation
			/// </summary>
			protected class Unsubsriber : IDisposable
			{
				/// <summary>
				/// List of observers of the simulation
				/// </summary>
				protected List<IObserver<SingleSimulation>> _observers;

				/// <summary>
				/// Observer of the simulation
				/// </summary>
				protected IObserver<SingleSimulation> _observer;

				/// <summary>
				/// Construct a new Unsubscriber
				/// </summary>
				/// <param name="observers">List of observers</param>
				/// <param name="observer">New observer</param>
				public Unsubsriber(List<IObserver<SingleSimulation>> observers, IObserver<SingleSimulation> observer)
				{
					this._observers = observers;
					this._observer = observer;
				}

				/// <summary>
				/// IDispoable interface implementer, Removes the observer from the list
				/// </summary>
				public void Dispose()
				{
					if (_observer != null && _observers.Contains(_observer))
						_observers.Remove(_observer);
				}
			}

			/// <summary>
			/// Ends the simulation and unsubscribes from it
			/// </summary>
			public void EndSimulation()
			{
				foreach (var observer in observers)
				{
					if (observers.Contains(observer))
						observer.OnCompleted();
				}

				observers.Clear();
			}
		}
	}
}
