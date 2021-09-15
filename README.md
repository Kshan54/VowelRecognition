# VowelRecognition
COMPILATION PROCEDURE:-
1)Open project testingvowel in folder 204101051_vowelRecognition and build it.
2).Ignore Conversion Warnings signed-unsigned / long int to int warnings while building.
3).Directly start Without degugging if build prompt is shown build it.

PROJECT DETAILS

1.2 projects one for testing and one for Reference File Generation.

  1a) testingvowel :- used for testing vowel test files / Live Testing.
 
  2a) reference_file :- used for generating average ci's of a vowel that one at a time and appending it to text file and create 25 * 12 

1.3 Runtime and Additional Details:

    1.3.1).Enter appropriate vowel as per menu options / Enter z to stop testing / Enter r for live testing.

    1.3.2).Code will print file name Vowel detected and Minimum distance From that reference vowel.

    1.3.3) If user wants to do testing for all 20 files per vowel testing as well refrence ones change 2 variable named total and count
       where Total denotes : total files per vowel to be tested:- 
       count denotes : Denotes start number of file so file number can be appended to string and tested continously.
       make Total = 20 and count =1 for testing reference as well as test files for testing that is 20 files per vowel.

    1.4.4).In Live testing Record Time is 4 seconds.

    1.4.5).DC shift is -0.0220778 calculated via separate silence File.

    1.4.6).For Steady Frames i am using 5 frames around "max ste" in Vowel (skipping starting 10 frames to avoid mic spike problem during live recording).  

    1.4.7).Accuracy For my testfiles  is 98.99 % .

    1.4.8).In Live testing A,I,E are showing 90% and above accuracy but error may occurs in O and U. 
