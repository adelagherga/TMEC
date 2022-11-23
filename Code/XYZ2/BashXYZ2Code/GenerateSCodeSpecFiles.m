/* 
To Run GenerateSCode: 
    1. Import ThesisCode Directory
    2. In parent directory of ThesisCode, enter 
            $ ABSOLUTE_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
       This determines the path to the parent directory containing ThesisCode
    3. Change current working directory to ThesisCode by entering 
            $ cd $ABSOLUTEPATH/ThesisCode
    4. In a terminal window, enter 
            $ source $ABSOLUTE_PATH/ThesisCode/GenerateSCode/BashGenerateSetsCode/GenerateSets.sh
    5. For desired m, in a terminal window, enter 
        $ GenerateSetsXYZ2 m

*/


{ GenerateSCode
    { GenerateSets.sh
        { GenerateSetsXYZ2.m 
            { nSets.m 
                { RecursiveSets.m } } } } 
            