import React, { useCallback, useMemo, useState } from 'react';
import { createRoot } from 'react-dom/client';
import { Editor } from 'ketcher-react';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import 'ketcher-react/dist/index.css';
import './style.css';

const exampleSmiles = '[OH-:1].[CH3:2][Br:3]>>[CH3:2][OH:1].[Br-:3]';

function App() {
  const structServiceProvider = useMemo(() => new StandaloneStructServiceProvider(), []);
  const [ketcher, setKetcher] = useState(null);
  const [status, setStatus] = useState('Ready');

  const loadExample = useCallback(async () => {
    if (!ketcher) {
      return;
    }
    await ketcher.setMolecule(exampleSmiles);
    setStatus('Loaded atom-mapped SN2 reaction');
  }, [ketcher]);

  const exportKet = useCallback(async () => {
    if (!ketcher) {
      return;
    }
    const ket = await ketcher.getKet();
    window.localStorage.setItem('codex-chem-last-ket', ket);
    setStatus(`Exported KET to localStorage (${ket.length} chars)`);
  }, [ketcher]);

  const exportSmiles = useCallback(async () => {
    if (!ketcher) {
      return;
    }
    const smiles = await ketcher.getSmiles();
    window.localStorage.setItem('codex-chem-last-smiles', smiles);
    setStatus(`Exported SMILES: ${smiles}`);
  }, [ketcher]);

  return (
    <div className="app-shell">
      <div className="toolbar">
        <button type="button" onClick={loadExample} disabled={!ketcher}>
          Load SN2
        </button>
        <button type="button" onClick={exportKet} disabled={!ketcher}>
          Save KET
        </button>
        <button type="button" onClick={exportSmiles} disabled={!ketcher}>
          Save SMILES
        </button>
        <span>{status}</span>
      </div>
      <div className="editor-frame">
        <Editor
          staticResourcesUrl="/"
          structServiceProvider={structServiceProvider}
          disableMacromoleculesEditor
          onInit={(instance) => {
            window.ketcher = instance;
            setKetcher(instance);
            setStatus('Ketcher initialized');
          }}
          errorHandler={(message) => setStatus(`Ketcher error: ${message}`)}
        />
      </div>
    </div>
  );
}

createRoot(document.getElementById('root')).render(<App />);
