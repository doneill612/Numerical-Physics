using System;
using System.Collections.Generic;

namespace Newtons_Method
{
    class ZeroFinder
    {
        //Constant Declarations
        static readonly double g = 9.81;
        static readonly double R = 8.32;
        static readonly double S = 1e-4;
        static readonly double P = 200;
        static readonly double T = 300;
        static readonly double l = 0.1;
        static readonly double x = P * S * l / (R * T);
        static readonly double m = 6e-3; 
        static readonly double omega_sqrd = g / l;
        static readonly double gamma = (2 * x * R) / (m * l * l);
        static readonly double theta0 = 1;
        
        /// <summary>
        /// A delegate that serves as a mathematical function of a single variable.
        /// </summary>
        /// <param name="x">function parameter</param>
        /// <returns>f(x)</returns>
        delegate double Delegate(double x);
        
        /// <summary>
        /// Contains a mapping of string keys to delegated functions.
        /// </summary>
        static Dictionary<string, Delegate> FunctionMap = new Dictionary<string, Delegate>();    

        // Application entry point
        static void Main(string[] args)
        {
            Init();

            double zero6_b = FindZero("sin", 0.1, 0.001);

            double zero6_c = FindZero("sin", 1.55, 0.001);

            Console.WriteLine("Response 6.b) " + 
                (zero6_b == double.MinValue ? "No root found" : zero6_b.ToString()));
            Console.WriteLine("Response 6.c) " +
                (zero6_c == double.MinValue ? "No root found" : zero6_c.ToString()));

            double zero6_d = FindZero("tanh", -4, 0.000001);

            Console.WriteLine("Response 6.d) " +
                (zero6_d == double.MinValue ? "No root found" : zero6_d.ToString()));

            double zero7 = FindZero("landau", 0.5, 0.00001);

            Console.WriteLine("Response 7) " +
                (zero7 == double.MinValue ? "No root found" : zero7.ToString()));

            Console.Read();

        }

        
        /// <summary>
        /// Initialize the FunctionMap with the functions in question.
        /// </summary>
        static void Init()
        {
            FunctionMap.Add("sin", Math.Sin);
            FunctionMap.Add("der_sin", Math.Cos);
            FunctionMap.Add("tanh", Tanh);
            FunctionMap.Add("der_tanh", TanhDeriv);
            FunctionMap.Add("landau", Landau);
            FunctionMap.Add("der_landau", LandauDeriv);
        }


        /// <summary>
        /// Uses Newton's Method to recursively find the zero of a particular function.
        /// </summary>
        /// <param name="functionKey">The key for the function mapping</param>
        /// <param name="xStart">The desired initial value for the interation</param>
        /// <param name="epsilon">The desired precision</param>
        /// <returns></returns>
        static double FindZero(string functionKey, double xStart, double epsilon)
        {
            Console.WriteLine("NEWTON'S METHOD: Finding zeros of function matching key: " + functionKey + "(x)");
            Console.WriteLine("Input parameters: x0 = " + xStart + ", eps = " + epsilon);

            // Apply "der_" marker to produce the key to the derivative function in the map
            string derivKey = string.Join("der_", new object[] { string.Empty, functionKey });

            // Fetch functions from the map
            Delegate f = FunctionMap[functionKey];
            Delegate df = FunctionMap[derivKey];

            // Define: xNext -> The next x value after the iteration.
            //         xCurrent -> The current x value in the iteration
            double xNext = 0;
            double xCurrent = xStart;

            // Define: zero -> The current zero location calculated after the iteration
            //         currentEpsilon -> The current level of precision after the iteration
            double zero = 0; 
            double currentEpsilon = xNext - xCurrent;

            // Keep track of the iteration, and define the maximum number of iterations allowed.
            int iteration = 0, maxIteration = 10000;

            // Begin iterations -> Continue until desired precision is bested, or iteration limit is reached.
            while (Math.Abs(currentEpsilon) > epsilon && iteration <= maxIteration)
            {
                xNext = xCurrent - ( f(xCurrent) / df(xCurrent) );

                currentEpsilon = xNext - xCurrent;
                zero = xNext;

                if (Math.Abs(xNext - xCurrent) > epsilon)
                {
                    xCurrent = xNext;
                }

                iteration++;
            }

            // Return the last found zero location.
            return iteration > maxIteration ? double.MinValue : zero; 
        }


        /// <summary>
        /// A function implementation for problem (6.d)
        /// </summary>
        /// <param name="x">The function parameter</param>
        /// <returns></returns>
        static double Tanh(double x)
        {
            return (1 / (2 - Math.Tanh(x - 1)));
        }

        /// <summary>
        /// The derivative of the function in problem (6.d) with respect to
        /// the function parameter.
        /// </summary>
        /// <param name="x">The function parameter</param>
        /// <returns></returns>
        static double TanhDeriv(double x)
        {
            double sech = 1 / Math.Cosh(1 - x);
            double sech2 = sech * sech;

            double tanh = Math.Tanh(1 - x) + 2;
            double tanh2 = tanh * tanh;

            double result = sech2 / tanh2;
            return result;
        }

        /// <summary>
        /// A Landau function implementation
        /// </summary>
        /// <param name="x">The function parameter</param>
        /// <returns></returns>
        static double Landau(double x)
        {
            return (
                ((Math.Sin(x)/x) * 
                (theta0*theta0 - x*x)) - 
                ((gamma * T)/(omega_sqrd))
            );
        }
        /// <summary>
        /// The derivative of the Landau function with respect to the function parameter
        /// </summary>
        /// <param name="x">The function parameter</param>
        static double LandauDeriv(double x)
        {
            return (
                (x * (theta0*theta0 - x*x) * Math.Cos(x) - 
                ((theta0*theta0 + x*x) * Math.Sin(x))) / (x*x)
            );
        }
    }
}
