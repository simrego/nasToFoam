/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------

Application
    nasToFoam

Group
    grpMeshConversionUtilities

Description
    nastran dat format mesh conversion.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "IFstream.H"
#include "StringStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Dat file format.
enum struct FORMAT : char
{
    FREE = 0,
    SMALL = 8,
    LARGE = 16
};
// File format, small by default.
FORMAT format = FORMAT::SMALL;

// Last extracted entry kw
string entryBuff;

// Buffer for last comment with it's line number. (patch, cellZone names)
word commentBuffer;
label commentLine;

// Get the next entry and store it in the buffer (fwd)
string& getEntry(IFstream& is);

// Ignore everything until we find the "BEGIN BULK" entry
bool findBulk(IFstream& is)
{
    string line;
    while (is.good())
    {
        is.getLine(line);
        if (line.starts_with("BEGIN BULK")) return true;
    }
    return false;
}

// Process a commented line.
// e.g. NX nastran write the property card name as comment before the entry.
// So now just store the last word from the comment which is the name.
void processCommentedLine(IFstream& is)
{
    string line;
    is.getLine(line);
    // Yep, this is not too safe...
    size_t lastSpace = line.find_last_of(' ');
    if (lastSpace < line.size())
    {
        commentBuffer = line.substr(lastSpace + 1);
        commentLine = is.lineNumber();
    }
}

// Skip everything on this line, and next lines if we have a multiline entry.
// The only problem if we get the '\n' from the stream. In this case the
// next line will be skipped which shouldn't be...
// WARNING: Make sure the '\n' is kept in the stream.
void finishEntry(IFstream& is)
{
    string line;
    while (true)
    {
        is.getLine(line);

        // Remove stupid carriage return char...
        line.erase(std::remove(line.begin(), line.end(), 13),
                 line.end());

        // If the end of the line is '+', the next line is the same entry.
        if (!line.ends_with('+')) break;
    }

}

// Extract the next column from the file.
string getColumn(IFstream& is, char size = char(format))
{
    if (size == 0) // Increase for free format to avoid overflow
    {
        size = 63;
    }
 
    label i = 0;
    char buff[63];  // More than enough, get size no chars

    while (i < size)
    {
        is.get(buff[i]);
        // Stop at new line, or ',' for free format.
        if (buff[i] == '\n')
        {
            i++;
            break;
        }
        else if (buff[i] == ',') // only possible in free format
        {
            // DOnt increase i, so we will change it to '\0'
            break;
        }
        i++;
    }
    buff[i] = '\0';

    // Continue on new line if end with "+\n" and next start with + or *
    if (buff[0] == '+' && buff[i-1] == '\n' &&
        (is.peek() == '+' || is.peek() == '*'))
    {
        if (format == FORMAT::FREE)
        {
            char c;
            while (is.get(c) && c != ',') {}
        }
        else is.readRaw(buff, 8);    // Ignore 1st column

        // Do it again.
        return getColumn(is, size);
    }

    // In free format '\n' is removed from the stream but we will need it
    // to finish the entry. Put it back IF we have it.
    if (format == FORMAT::FREE && buff[i-1] == '\n')
    {
        buff[i-1] = '\0';
        is.putback('\n');
    }

    string data = string(buff);
    // Remove whitespaces
    data.erase(std::remove(data.begin(), data.end(), ' '),
                 data.end());

    // Remove stupid carriage return char...
    data.erase(std::remove(data.begin(), data.end(), 13),
             data.end());

    return data;
}

string& getEntry(IFstream& is)
{
    // Skip comments
    while (is.good() && is.peek() == '$') processCommentedLine(is);

    switch (format)
    {
    case FORMAT::SMALL:
    case FORMAT::LARGE:
        entryBuff = getColumn(is, 8);
        break;
    case FORMAT::FREE:
        entryBuff = getColumn(is, 0);
        break;
    }

    if (entryBuff.ends_with('*')) // We have multiline, ignore *.
        entryBuff.removeEnd('*');

    return entryBuff;
}

label getLabel(IFstream& is)
{
    return readLabel(getColumn(is));
}

// Scientific notation sucks in nastran...
// Sometimes we have E, sometimes we don't... ?!?!
scalar getScalar(IFstream& is)
{
    string data = getColumn(is);
    // This is stupid. And probably costly.
    label id = data.find('+', 1); // contains e+
    if (id == -1) id = data.find('-', 1); // contains e-
    if (id != -1 && (data[id-1] != 'e' && data[id-1] != 'E'))
    {
        data = data.substr(0, id) + "e" + data.substr(id);
    }

    return readFloat(data);
}

// Read points until we find some different entry.
// GRID card format, where CP is ignored:
// GRID   ID   CP   X  Y  Z  ...
void readPoints
(
    IFstream& is,
    DynamicList<point>& points,
    DynamicList<label>& pointIDs
)
{
    Info<< "\tReading \"GRID\" entries." << endl;

    label pointi;
    while (true)
    {
        pointi = getLabel(is);
        if (pointi >= pointIDs.size())
        {
            pointIDs.resize(pointi + 1, -1);
        }
        pointIDs[pointi] = points.size();

        // Ignore CP column...
        getColumn(is);
        // Get the 3 coordinate
        point pt;
        pt[0] = getScalar(is);
        pt[1] = getScalar(is);
        pt[2] = getScalar(is);
        points.append(pt);
        // Ignore the other entries
        finishEntry(is); //is.getLine(nullptr);

        // No more GRID entry.
        if (!getEntry(is).starts_with("GRID")) break;
    }
}

// Read cells of a given type.
template<cellModel::modelType TYPE>
void readCell
(
    string name,
    IFstream& is,
    DynamicList<cellShape>& cellVerts,
    DynamicList<DynamicList<label>>& cellPropIDs,
    DynamicList<label>& pointIDs
)
{
    Info<< "\tReading " << name << " entries." << endl;
    
    const cellModel& model = cellModel::ref(TYPE);
    label cellPropID;
    while (entryBuff == name)
    {
        getLabel(is);   // cell ID
        cellPropID = getLabel(is);
        if (cellPropID >= cellPropIDs.size())
        {
            cellPropIDs.resize(cellPropID + 1);
        }
        cellPropIDs[cellPropID].append(cellVerts.size());
        
        labelList verts(model.nPoints());
        forAll(verts, i)
        {
            verts[i] = pointIDs[getLabel(is)];
        }
        cellVerts.append(cellShape(model, verts, true));

        finishEntry(is);
        getEntry(is);
    }
}

template<label nPoints>
void readFaces
(
    string name,
    IFstream& is,
    DynamicList<DynamicList<face>>& patches,
    DynamicList<label>& pointIDs
)
{
    label patchi;
    while (entryBuff == name)
    {
        getLabel(is); // ignore ID
        patchi = getLabel(is);
        if (patchi >= patches.size())
        {
            patches.resize(patchi + 1);
        }

        face fVerts(nPoints);
        forAll(fVerts, i)
        {
            fVerts[i] = pointIDs[getLabel(is)];
        }
        patches[patchi].append(fVerts);

        finishEntry(is); // Get the end of the line
        getEntry(is);
    }
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert nastran dat file to OpenFOAM. Only for small format and Meters."
    );
    argList::noParallel();
    argList::addArgument(".dat file");
    argList::addOption(
        "format",
        "word",
        "Input format. Options: small(default), large, free"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    if (args.found("format"))
    {
        word inpForm = args.getOrDefault<word>("format", "small");
        if (inpForm == "small") 
        {
            // Do nothing. Initialised to small.
        }
        else if (inpForm == "large")
        {
            format = FORMAT::LARGE;
        }
        else if (inpForm == "free")
        {
            format = FORMAT::FREE;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown format: " << inpForm << "."
                << exit(FatalError);
        }
    }

    const auto datName = args.get<fileName>(1);
    IFstream inFile(datName);

    if (!inFile.good())
    {
        FatalErrorInFunction
            << "Cannot open file " << datName
            << exit(FatalError);
    }

    if (!findBulk(inFile))
    {
        FatalErrorInFunction
            << "Cannot find \"BEGIN BULK\" entry."
            << exit(FatalError);
    }

    // Points
    DynamicList<point> points;
    // Nastran indexing
    DynamicList<label> pointIDs;
    // Cell vertices
    DynamicList<cellShape> cellVerts;
    // Cell property IDs
    DynamicList<DynamicList<label>> cellPropIDs;
    // Patch IDs
    DynamicList<DynamicList<face>> patches;
    DynamicList<label> patchIDs;
    // Porperty card names
    Map<word> cellZonePropNames;
    Map<word> patchPropNames;

    // Read the next entry into the buffer.
    getEntry(inFile);

    Info<< "Start reading file." << endl;
    while (inFile.good())
    {
        if (entryBuff == "GRID")
        {
            readPoints(inFile, points, pointIDs);
        }
        else if (entryBuff == "CTETRA")
        {
            readCell<cellModel::modelType::TET>(
                "CTETRA", inFile, cellVerts, cellPropIDs, pointIDs
            );
        }
        else if (entryBuff == "CPYRAM")
        {
            readCell<cellModel::modelType::PYR>(
                "CPYRAM", inFile, cellVerts, cellPropIDs, pointIDs
            );
        }
        else if (entryBuff == "CHEXA")
        {
            readCell<cellModel::modelType::HEX>(
                "CHEXA", inFile, cellVerts, cellPropIDs, pointIDs
            );
        }
        else if (entryBuff == "CTRIA3")
        {
            readFaces<3>(
                "CTRIA3", inFile, patches, pointIDs
            );
        }
        else if (entryBuff == "PSOLID")
        {
            // Cell Zone names
            label czI = getLabel(inFile);

            if (commentLine == inFile.lineNumber())
            {
                cellZonePropNames.insert(czI, commentBuffer);
            }
            else
            {
                cellZonePropNames.insert(czI, "");
            }
            finishEntry(inFile);
            getEntry(inFile);
        }
        else if (entryBuff == "PSHELL")
        {
            // Patches
            label patchI = getLabel(inFile);
            if (commentLine == inFile.lineNumber())
            {
                patchPropNames.insert(patchI, commentBuffer);
            }
            else
            {
                patchPropNames.insert(patchI, "");
            }
            finishEntry(inFile);
            getEntry(inFile);
        }
        else if (entryBuff == "ENDDATA")
        {
            Info<< "Finished reading file." << endl;
            break;
        }
        else
        {
            FatalErrorInFunction
                << "Cannot process keyword: \"" << entryBuff.c_str()
                << "\", on line " << inFile.lineNumber() << "."
                << exit(FatalError);
        }
    }

    DynamicList<faceList> patchFaces;
    wordList patchNames;
    label unnamedPatchN = 0;
    forAll(patchPropNames, i)
    {
        label propI = patchPropNames.toc()[i];
        word propName = patchPropNames.at(propI);
        if (patches[propI].size())
        {
            patchFaces.append(std::move(patches[propI]));
            if (propName.empty())
            {
                patchNames.append("patch_" + std::to_string(unnamedPatchN++));
            }
            else
            {
                patchNames.append(propName);
            }
        }
    }

    Info<< "Constructing the mesh." << endl;
    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        pointField(points),
        cellVerts,
        patchFaces,
        patchNames,
        wordList(patchNames.size(), polyPatch::typeName),
        "defaultFaces",
        polyPatch::typeName,
        wordList()
    );

    
    if (cellPropIDs.size())
    {
        Info<< "Adding cell zones." << endl;

        List<cellZone*> cZones;
        label unnamedCellZoneN = 0;
        forAll(cellZonePropNames, i)
        {
            label propI = cellZonePropNames.toc()[i];
            word propName = cellZonePropNames.at(propI);
            if (propName.empty())
                propName = "cellZone_" + std::to_string(unnamedCellZoneN++);
            cZones.append
            (
                new cellZone
                (
                    propName,
                    cellPropIDs[propI],
                    i,
                    mesh.cellZones()
                )
            );
        }

        mesh.addZones(List<pointZone*>(), List<faceZone*>(), cZones);
    }

    Info<< endl;
    Info<< "Mesh information:" << endl
        << "Number of points: " << mesh.nPoints() << endl
        << "Number of faces: " << mesh.nFaces() << endl
        << "Number of cells: " << mesh.nCells() << endl
        << "Found patch names:" << endl;
    forAll(patchNames, i) Info<< "\t" << patchNames[i] << endl;
    if (mesh.cellZones().size())
    {
        Info<< "Found cell zones:" << endl;
        forAll(mesh.cellZones(), i)
        {
            Info<< "\t" << mesh.cellZones()[i].name() << endl;
        }
    }
    Info<< endl;

    mesh.removeFiles();
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
