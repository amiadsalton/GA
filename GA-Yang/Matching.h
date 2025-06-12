
using System;
using System.IO;

namespace EvolutionaryComputing
{
    //******************************************************************************
    // Name: CMatching
    // Role: Reperesents matching preferences
    // Responsiblities: 
    //      1. Knows to retrieve and save initial preferences from file
    //      2. Knows to create new (randomized) preferences set
    //      3. Knows to evalute a group of (women, man and dog). The evaluation values is 
    //         the average distance between selection to preference.
    //******************************************************************************
    public class CMatching
    {
        // Enumaration and constants:
        public enum EIndividualType { eWoman = 0, eMan, eDog };
        public const uint DATA_SIZE = 50;

        //====================
        // Method: CMatching
        // Responsibility: Constructor
        //====================
        public CMatching()
        {
            // Initiates preferences list
            mPreferences = new CIndividualPreferences[3, DATA_SIZE];
            for (uint type = (uint)EIndividualType.eWoman; type < 3; type++)
            {
                // For each individual of that type
                for (uint individual = 0; individual < DATA_SIZE; individual++)
                {
                    mPreferences[type, individual] = new CIndividualPreferences();
                }
            }
        }

        //====================
        // Method: Init
        // Responsibility: Init preferences tabel
        // Method: Read from file, and if file doesn't exist create new randomized list
        //====================
        public void Init(bool init)
        {
            if (init == false)
            {
                try
                {
                    // Try to read from file
                    StreamReader reader = new StreamReader("Preferences.txt");
                    FromStream(reader);
                    reader.Close();
                    return;
                }
                catch
                {
                }
            }

            // File doesn't exists (or need ot build a new list)
            BuildPreferencesData();     // Build the list

            // Save it into preferences file
            StreamWriter writer = new StreamWriter("Preferences.txt");
            writer.Write(ToString());
            writer.Close();
        }

        //====================
        // Method: Evaluate
        // Responsibility: Evaluate a match (women, men, dog)
        // Method: average distance between selection to preference
        //====================
        public uint Evaluate(uint[] iGroup)
        {
            uint sum = 0;
            for (uint i = 0; i < 3; i++)
            {
                CIndividualPreferences pref = mPreferences[i, iGroup[i]];
                sum += pref.mFirstPreferences[iGroup[(i + 1) % 3]] + pref.mSecondPreferences[iGroup[(i + 2) % 3]];
            }
            return (uint)((float)sum / (3 * (3 - 1)) + 0.5);
        }

        //====================
        // Method: BuildPreferencesData
        // Responsibility: Builds a new preferences table
        // Method: Build a pool of all items and randomly select item at each step
        //====================
        private void BuildPreferencesData()
        {
            Random random = new Random(unchecked((int)DateTime.Now.Ticks));
            // For each type
            for (uint type = (uint)EIndividualType.eWoman; type < 3; type++)
            {
                // For each individual of that type
                for (uint individual = 0; individual < DATA_SIZE; individual++)
                {
                    // Define a pool and allocate for each element preference
                    uint[] firstPool = new uint[DATA_SIZE];
                    uint[] secondPool = new uint[DATA_SIZE];
                    // Fill pool
                    for (uint i = 0; i < DATA_SIZE; i++)
                    {
                        firstPool[i] = secondPool[i] = i;
                    }

                    // Allocate preferences from pool for each index
                    for (uint i = 0; i < DATA_SIZE; i++)
                    {
                        uint firstPoolIndex = (uint)random.Next(0, (int)(DATA_SIZE - i));
                        uint secondPoolIndex = (uint)random.Next(0, (int)(DATA_SIZE - i));
                        mPreferences[type, individual].mFirstPreferences[i] = firstPool[firstPoolIndex];
                        mPreferences[type, individual].mSecondPreferences[i] = secondPool[secondPoolIndex];

                        // Delete element from pool
                        for (uint j = firstPoolIndex; j < DATA_SIZE - i - 1; j++)
                        {
                            firstPool[j] = firstPool[j + 1];
                        }
                        for (uint j = secondPoolIndex; j < DATA_SIZE - i - 1 ; j++)
                        {
                            secondPool[j] = secondPool[j + 1];
                        }
                    }
                }
            }
        }

        //====================
        // Method: ToString
        // Responsibility: create a string out of all preferences
        //====================
        public override string ToString()
        {
            string result=new string(' ',0);
            for (uint type = (uint)EIndividualType.eWoman; type < 3; type++)
            {
                // For each individual of that type
                for (uint individual = 0; individual < DATA_SIZE; individual++)
                {
                    result += mPreferences[type, individual].ToString();
                    result += "\n\n";
                }
                result += "\n\n";
            }
            return result;
        }

        //====================
        // Method: FromStream
        // Responsibility: Read preferences from file
        //====================
        private void FromStream(StreamReader input)
        {
            for (uint type = (uint)EIndividualType.eWoman; type < 3; type++)
            {
                // For each individual of that type
                for (uint individual = 0; individual < DATA_SIZE; individual++)
                {
                    mPreferences[type, individual].FromStream(input);
                }
                string line = input.ReadLine();
                line = input.ReadLine();
            }
        }

        //====================
        // Class: CIndividualPreferences
        // Role: Internal class which represents individual preferences
        //====================
        private class CIndividualPreferences
        {
            // constructor
            public CIndividualPreferences()
            {
                mFirstPreferences = new uint[DATA_SIZE];
                mSecondPreferences = new uint[DATA_SIZE];
            }

            // Create string out of this preferences
            public override string ToString()
            {
                string result = new string(' ',0);
                for (uint i = 0; i < DATA_SIZE; i++)
                {
                    result += mFirstPreferences[i].ToString();
                    result += " ";
                }
                result += "\n";
                for (uint i = 0; i < DATA_SIZE; i++)
                {
                    result += mSecondPreferences[i].ToString();
                    result += " ";
                }
                result += "\n";
                return result;
            }

            // Read indivisuals preferences out of file
            public void FromStream(StreamReader input)
            {
                string line=input.ReadLine();
                string[] preferences = line.Split(' '); 
                for (uint i = 0; i < DATA_SIZE; i++)
                {
                    mFirstPreferences[i] = uint.Parse(preferences[i]);
                }
                line = input.ReadLine();
                preferences = line.Split(' ');
                for (uint i = 0; i < DATA_SIZE; i++)
                {
                    mSecondPreferences[i] = uint.Parse(preferences[i]);
                }
                line = input.ReadLine();
                line = input.ReadLine();
            }

            // first type preferences (e.g. men for woman)
            public uint[] mFirstPreferences;
            // second type preferences (e.g. dogs for woman)
            public uint[] mSecondPreferences;
        };

        // All preferences
        private CIndividualPreferences[,]   mPreferences;

        // Maximum possible evaluation value
        public static uint sMaxEvaluation = DATA_SIZE;
    }
}