import React, { useCallback, useEffect, useMemo, useState } from 'react';
import { createRoot } from 'react-dom/client';
import './raphael-shim.js';
import { Editor } from 'ketcher-react';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import 'ketcher-react/dist/index.css';
import './style.css';

const exampleSmiles = '[OH-:1].[CH3:2][Br:3]>>[CH3:2][OH:1].[Br-:3]';
const params = new URLSearchParams(window.location.search);
const reviewConfig = {
  sessionId: params.get('session') || params.get('session_id'),
  token: params.get('token'),
  api: (params.get('api') || window.location.origin).replace(/\/$/, ''),
};
const hasReviewSession = Boolean(reviewConfig.sessionId && reviewConfig.token && reviewConfig.api);
const terminalStatuses = new Set(['confirmed', 'modified', 'skipped', 'error']);

class EditorErrorBoundary extends React.Component {
  constructor(props) {
    super(props);
    this.state = { error: null };
  }

  static getDerivedStateFromError(error) {
    return { error };
  }

  componentDidCatch(error, info) {
    window.codexKetcherError = {
      message: error?.message || String(error),
      stack: error?.stack || '',
      componentStack: info?.componentStack || '',
    };
  }

  render() {
    if (this.state.error) {
      const message = this.state.error?.message || String(this.state.error);
      return (
        <div className="editor-error" role="alert">
          <strong>Ketcher editor failed to render.</strong>
          <pre>{message}</pre>
        </div>
      );
    }
    return this.props.children;
  }
}

function itemTitle(item, index) {
  if (!item) {
    return '';
  }
  return item.label || item.id || `Item ${index + 1}`;
}

function App() {
  const structServiceProvider = useMemo(() => new StandaloneStructServiceProvider(), []);
  const [ketcher, setKetcher] = useState(null);
  const [status, setStatus] = useState('Ready');
  const [session, setSession] = useState(null);
  const [currentIndex, setCurrentIndex] = useState(0);
  const [posting, setPosting] = useState(false);

  const apiRequest = useCallback(async (path, options = {}) => {
    const separator = path.includes('?') ? '&' : '?';
    const response = await fetch(
      `${reviewConfig.api}${path}${separator}token=${encodeURIComponent(reviewConfig.token)}`,
      {
        ...options,
        headers: {
          'Content-Type': 'application/json',
          'X-Review-Token': reviewConfig.token,
          ...(options.headers || {}),
        },
      },
    );
    const payload = await response.json();
    if (!response.ok) {
      const message = payload.warnings?.join('; ') || `HTTP ${response.status}`;
      throw new Error(message);
    }
    return payload;
  }, []);

  const refreshSession = useCallback(async () => {
    if (!hasReviewSession) {
      return;
    }
    const payload = await apiRequest(`/api/review-sessions/${encodeURIComponent(reviewConfig.sessionId)}`);
    setSession(payload);
    const pendingIndex = payload.items.findIndex((item) => item.status === 'pending');
    if (pendingIndex >= 0) {
      setCurrentIndex(pendingIndex);
    }
    setStatus(payload.status === 'pending' ? 'Review session loaded' : `Session ${payload.status}`);
  }, [apiRequest]);

  useEffect(() => {
    if (!hasReviewSession) {
      return;
    }
    refreshSession().catch((error) => setStatus(`Review load error: ${error.message}`));
  }, [refreshSession]);

  const items = session?.items || [];
  const currentItem = items[currentIndex] || null;
  const completedCount = items.filter((item) => terminalStatuses.has(item.status)).length;

  useEffect(() => {
    if (!hasReviewSession || !ketcher || !currentItem || currentItem.status !== 'pending') {
      return;
    }
    let cancelled = false;
    ketcher
      .setMolecule(currentItem.input_smiles || currentItem.smiles)
      .then(() => {
        if (!cancelled) {
          setStatus(`Loaded ${itemTitle(currentItem, currentIndex)}`);
        }
      })
      .catch((error) => {
        if (!cancelled) {
          setStatus(`Ketcher load error: ${error.message}`);
        }
      });
    return () => {
      cancelled = true;
    };
  }, [currentIndex, currentItem, ketcher]);

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

  const submitCurrent = useCallback(
    async (submissionStatus) => {
      if (!currentItem || posting || (submissionStatus !== 'skipped' && !ketcher)) {
        return;
      }
      setPosting(true);
      try {
        const body = { status: submissionStatus };
        if (submissionStatus !== 'skipped') {
          try {
            body.reviewed_smiles = await ketcher.getSmiles();
            body.molfile = await ketcher.getMolfile();
            const ket = await ketcher.getKet();
            window.localStorage.setItem(`codex-chem-review-ket-${currentItem.id}`, ket);
          } catch (error) {
            body.status = 'error';
            body.reviewed_smiles = currentItem.input_smiles || currentItem.smiles;
            body.warnings = [`Ketcher export failed: ${error.message}`];
          }
        }
        const payload = await apiRequest(
          `/api/review-sessions/${encodeURIComponent(reviewConfig.sessionId)}/items/${encodeURIComponent(currentItem.id)}`,
          {
            method: 'POST',
            body: JSON.stringify(body),
          },
        );
        setSession(payload);
        const nextPending = payload.items.findIndex((item) => item.status === 'pending');
        if (nextPending >= 0) {
          setCurrentIndex(nextPending);
          setStatus(`Saved ${itemTitle(currentItem, currentIndex)}`);
        } else {
          setStatus(`Session ${payload.status}`);
        }
      } catch (error) {
        setStatus(`Submit error: ${error.message}`);
      } finally {
        setPosting(false);
      }
    },
    [apiRequest, currentIndex, currentItem, ketcher, posting],
  );

  const completeSession = useCallback(async () => {
    if (posting) {
      return;
    }
    setPosting(true);
    try {
      const payload = await apiRequest(
        `/api/review-sessions/${encodeURIComponent(reviewConfig.sessionId)}/complete`,
        { method: 'POST', body: '{}' },
      );
      setSession(payload);
      setStatus(`Session ${payload.status}`);
    } catch (error) {
      setStatus(`Complete error: ${error.message}`);
    } finally {
      setPosting(false);
    }
  }, [apiRequest, posting]);

  const cancelSession = useCallback(async () => {
    if (posting) {
      return;
    }
    setPosting(true);
    try {
      const payload = await apiRequest(
        `/api/review-sessions/${encodeURIComponent(reviewConfig.sessionId)}/cancel`,
        { method: 'POST', body: '{}' },
      );
      setSession(payload);
      setStatus(`Session ${payload.status}`);
    } catch (error) {
      setStatus(`Cancel error: ${error.message}`);
    } finally {
      setPosting(false);
    }
  }, [apiRequest, posting]);

  const reviewToolbar = (
    <div className="toolbar">
      <div className="queue-summary">
        <strong>{session?.status || 'Loading'}</strong>
        <span>
          {completedCount}/{items.length}
        </span>
      </div>
      <button
        type="button"
        onClick={() => submitCurrent('confirmed')}
        disabled={!ketcher || !currentItem || currentItem.status !== 'pending' || posting || session?.status !== 'pending'}
      >
        Confirm
      </button>
      <button
        type="button"
        onClick={() => submitCurrent('modified')}
        disabled={!ketcher || !currentItem || currentItem.status !== 'pending' || posting || session?.status !== 'pending'}
      >
        Save Edit
      </button>
      <button
        type="button"
        onClick={() => submitCurrent('skipped')}
        disabled={!currentItem || currentItem.status !== 'pending' || posting || session?.status !== 'pending'}
      >
        Skip
      </button>
      <button type="button" onClick={completeSession} disabled={!session || posting || session?.status !== 'pending'}>
        Finish
      </button>
      <button type="button" onClick={cancelSession} disabled={!session || posting || session?.status !== 'pending'}>
        Cancel
      </button>
      <span>{status}</span>
    </div>
  );

  const demoToolbar = (
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
  );

  return (
    <div className="app-shell">
      {hasReviewSession ? reviewToolbar : demoToolbar}
      {hasReviewSession && (
        <div className="queue-panel">
          {items.map((item, index) => (
            <button
              className={index === currentIndex ? 'queue-item active' : 'queue-item'}
              disabled={posting || session?.status !== 'pending'}
              key={item.id}
              onClick={() => setCurrentIndex(index)}
              type="button"
            >
              <span>{itemTitle(item, index)}</span>
              <small>{item.status}</small>
            </button>
          ))}
        </div>
      )}
      <div className="editor-frame">
        <EditorErrorBoundary>
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
        </EditorErrorBoundary>
      </div>
    </div>
  );
}

createRoot(document.getElementById('root')).render(<App />);
