//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// File: G4ElementSelector
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
//
// Creation date: 2 April 2000
//
// Modifications:
// 18/08/2000  V.Ivanchenko Update description
// 17/05/2006  V.Ivanchenko Cleanup
// 02/10/2007  V.Ivanchenko Fixed typo in computation of Lambda-factor
//                          proposed by Victor Pec
//
//---------------------------------------------------------------------

#include "G4ElementSelector.hh"

#include "G4Material.hh"
#include "G4MuonMinus.hh"
#include "G4Nucleus.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

G4bool G4ElementSelector::fMuonicHydrogenTransfer = false;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ElementSelector::G4ElementSelector() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ElementSelector::~G4ElementSelector() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Element* G4ElementSelector::SelectZandA(const G4Track& track, G4Nucleus* target)
{
  // Fermi-Teller Z-low of mu- capture and exceptions
  // for halogens and oxigen.
  // N.C.Mukhopadhyay Phys. Rep. 30 (1977) 1.

  std::size_t i = 0;
  const G4Material* mat = track.GetMaterial();
  std::size_t numberOfElements = mat->GetNumberOfElements();
  const G4ElementVector* theElementVector = mat->GetElementVector();

  if (1 < numberOfElements)
  {
    if (numberOfElements > prob.size())
    {
      prob.resize(numberOfElements, 0.0);
    }

    const G4double* theAtomNumDensity = mat->GetAtomicNumDensityVector();

    G4double sum = 0.0;
    for (i = 0; i < numberOfElements; ++i)
    {
      G4int Z = (*theElementVector)[i]->GetZasInt();

      // Halogens
      if ((9 == Z) || (17 == Z) || (35 == Z) || (53 == Z) || (85 == Z))
      {
        sum += 0.66 * Z * theAtomNumDensity[i];

        // Oxigen
      }
      else if (8 == Z)
      {
        sum += 0.56 * Z * theAtomNumDensity[i];

        // Others
      }
      else
      {
        sum += Z * theAtomNumDensity[i];
      }
      prob[i] = sum;
    }

    sum *= G4UniformRand();
    for (i = 0; i < numberOfElements; ++i)
    {
      if (sum <= prob[i])
      {
        break;
      }
    }

    // Optional (off by default): transfer of the muon from hydrogen.
    // A mu-p atom is small and neutral; in a solid or a liquid it gives
    // its muon to a heavier nucleus at a rate of 1e10-1e11 /s, far above
    // the decay rate 4.55e5 /s (D.F.Measday, Phys. Rep. 354 (2001) 243).
    // A mu- selected on hydrogen is re-assigned to a Z>1 element drawn
    // with the weights above. The random number already drawn is
    // re-scaled, so no extra random number is used.
    if (fMuonicHydrogenTransfer && 1 == (*theElementVector)[i]->GetZasInt()
        && track.GetDefinition() == G4MuonMinus::MuonMinus())
    {
      const G4State state = mat->GetState();
      if (kStateSolid == state || kStateLiquid == state
          || (kStateUndefined == state && mat->GetDensity() >= 0.1 * CLHEP::g / CLHEP::cm3))
      {
        G4double wheavy = 0.0;
        for (std::size_t j = 0; j < numberOfElements; ++j)
        {
          if (1 < (*theElementVector)[j]->GetZasInt())
          {
            wheavy += prob[j] - ((0 < j) ? prob[j - 1] : 0.0);
          }
        }
        if (0.0 < wheavy)
        {
          G4double low = (0 < i) ? prob[i - 1] : 0.0;
          G4double x = wheavy * (sum - low) / (prob[i] - low);
          G4double acc = 0.0;
          for (std::size_t j = 0; j < numberOfElements; ++j)
          {
            if (1 < (*theElementVector)[j]->GetZasInt())
            {
              acc += prob[j] - ((0 < j) ? prob[j - 1] : 0.0);
              i = j;
              if (x <= acc)
              {
                break;
              }
            }
          }
        }
      }
    }
  }

  const G4Element* elm = (*theElementVector)[i];
  G4int Z = elm->GetZasInt();

  // select isotope
  const G4IsotopeVector* isv = elm->GetIsotopeVector();
  std::size_t ni = isv->size();
  i = 0;

  if (1 < ni)
  {
    const G4double* ab = elm->GetRelativeAbundanceVector();
    G4double y = G4UniformRand();
    for (i = 0; i < ni; ++i)
    {
      y -= ab[i];
      if (y <= 0.0)
      {
        break;
      }
    }
  }

  G4int A = elm->GetIsotope((G4int)i)->GetN();
  target->SetParameters(A, Z);

  return elm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ElementSelector::SetMuonicHydrogenTransfer(G4bool val)
{
  fMuonicHydrogenTransfer = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4ElementSelector::GetMuonicHydrogenTransfer()
{
  return fMuonicHydrogenTransfer;
}
