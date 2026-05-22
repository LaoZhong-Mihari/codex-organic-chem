import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';

function ketcherRaphaelRequireShim() {
  return {
    name: 'ketcher-raphael-require-shim',
    enforce: 'pre',
    transform(code) {
      if (!code.includes('require("raphael")') && !code.includes("require('raphael')")) {
        return null;
      }
      return {
        code: code
          .replaceAll('require("raphael")', 'globalThis.Raphael')
          .replaceAll("require('raphael')", 'globalThis.Raphael'),
        map: null,
      };
    },
  };
}

export default defineConfig({
  plugins: [ketcherRaphaelRequireShim(), react()],
});
